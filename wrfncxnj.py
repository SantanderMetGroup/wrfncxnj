#!/usr/bin/env python

import netCDF4 as ncdf
import logging
import string
import sys
from glob import glob
from options import get_options
import wrfncxnj_fun
from wrfncxnj_base import (WrfNcFiles, OfileHandler, Projection,
                           read_variable_table, WrfNcTime, strptime,
                           create_oncs, ParseTransform, get_oncvar)
log = logging.getLogger(__name__)


class DeprecationError(RuntimeError):
    pass


def check_options_consistency(opt):
    if opt.refdate:
        log.error("Error: This option has been deprecated, please use -u "
                  "--time-units instead")
        raise DeprecationError

    if not opt.OFILE and not opt.OUTPUT_PATTERN:
        raise IOError("Missing output file or pattern!")
    if opt.OFILE and opt.OUTPUT_PATTERN:
        raise IOError(
            "-o and --output-pattern cannot be used at the same time!")
    if opt.splitvars and not opt.OUTPUT_PATTERN:
        raise AssertionError(
            "Split variables has been activated, but an output "
            "pattern has not been specified. Use --output-pattern="
            " to specify one"
        )

    if not opt.splitvars and opt.OUTPUT_PATTERN:
        raise AssertionError("An output pattern is needed to use "
                             "--split-in-variables. Please set it with the "
                             "--output-pattern option")
    if not opt.splitvars and opt.selected_plevs:
        raise AssertionError("--split-in-vars and -sel-plevs must be used at "
                             "the same time")


class ExtractAndJoin(object):
    def __init__(self, opt, args):
        self.opt = opt
        self.args = args
        self.wrf_files_iterator = WrfNcFiles()
        if opt.geofile:
            self.wrf_files_iterator.geo = ncdf.Dataset(opt.geofile, "r")
        if opt.fullfile:
            self.wrf_files_iterator.full = ncdf.Dataset(opt.fullfile, "r")
        check_options_consistency(opt)
        self.ofh = OfileHandler(opt)
        self.ofh.set_tempdir()
        self.is_geofile = False
        self.files = self.get_input_files()
        self.wrf_files_iterator.loadFiles(self.files)
        self.requested_variable_names = opt.requested_variables.split(',')
        if not opt.vtable:
            self.opt.vtable = sys.argv[0].replace(".py", ".table")

        self.projection = self.get_projection()
        self.variables = read_variable_table(
            self.requested_variable_names,
            opt.vtable,
            self.projection,
            self.files[0],
            opt
        )
        refdate = opt.time_units.split()[2]
        self.wrfnctime = WrfNcTime(strptime(refdate), opt.time_units)

    def run(self):
        opt = self.opt
        #
        # Create output netCDF files
        #
        oncs, oncnames = create_oncs(
            self.variables,
            self.ofh,
            self.files,
            self.wrfnctime,
            self.projection,
            opt
        )
        self.ofh.set_oncnames(oncnames)
        #
        # Get coordinate variable objects from the output netCDF files.
        #
        if self.opt.discard:
            raise NotImplementedError
            # files = discard_suspect_files(files, opt.discard)
        #
        # Loop over files extracting variables and times
        #
        for wrfnc in self.wrf_files_iterator:
            oncvar = self.process_file(wrfnc, oncs)
        oncvar.sync()
        oncvar.close()
        if opt.ftimes:
            time1, time2 = opt.ftimes.split(",")
            self.ofh.filter_times(time1, time2)
        if opt.TEMPDIR:
            self.ofh.move_output()
        self.wrf_files_iterator.closeCommon()
        print "Extract and join finished successfully"

    def process_file(self, wrfnc, oncs):
        opt = self.opt
        log.debug(wrfnc)
        if not opt.singlerec:
            self.wrfnctime.checkStep(wrfnc)
        times = self.wrfnctime.getTimes(wrfnc, opt.singlerec)
        log.debug(self.wrfnctime)
        for varname in self.variables:
            varobj = self.variables[varname]
            oncvar = self.process_variable_in_file(varname, varobj, wrfnc, oncs,
                                                   times)
        self.wrfnctime.cycle()
        return oncvar

    def process_variable_in_file(self, varname, varobj, wrfnc, oncs, times):
        opt = self.opt
        wnt = self.wrfnctime
        if not opt.splitvars:
            onc = oncs
        else:
            onc = oncs[varname]
        if not opt.quiet:
            print "Processing var %s" % varname

        function_name = "compute_%s" % varname
        if hasattr(wrfncxnj_fun, function_name):
            log.debug("Using the %s function" % function_name)
            process_func = getattr(wrfncxnj_fun, function_name)
            oncvar, copyval = process_func(varobj, onc, wrfnc, wnt)
            oncvar.set_ncvars(varobj)
            oncvar.assign_values(
                copyval[:wnt.nrec].astype('f'),
                wnt.iini,
                wnt.iend,
                times
            )
            oncvar.sync()
        elif varobj.transform:
            parse_transform = ParseTransform(varobj.transform, self.opt)
            log.debug(parse_transform)
            vardic = {}
            for var in parse_transform.variables:
                vardic[var] = wrfnc.current.variables[var][:]
            oncvar, copyval = parse_transform.execute(
                varobj,
                onc,
                wrfnc,
                wnt,
                vardic
            )
            oncvar.set_ncvars(varobj)
            oncvar.assign_values(copyval[:wnt.nrec].astype('f'), wnt.iini,
                                 wnt.iend, times)
            oncvar.sync()
        else:
            incvar = wrfnc.current.variables[varname]
            oncvar = get_oncvar(varobj, incvar, onc)
            oncvar.set_ncvars(varobj)
            oncvar.assign_values(incvar[:wnt.nrec].astype('f'), wnt.iini,
                                 wnt.iend, times)
            oncvar.sync()
        return oncvar

    def get_input_files(self):
        opt = self.opt
        if opt.globfiles:
            files = glob(opt.globfiles)
            files.sort()
        elif opt.filelist:
            files = map(string.strip, open(opt.filelist, "r").readlines())
        else:
            files = self.args
        #
        # Look for geofile
        #
        if not files and self.wrf_files_iterator.geo:
            print "No input files provided."
            print "Trying to find the variables in the geo_em file provided."
            files = [opt.geofile, ]
            self.is_geofile = True
        return files

    def get_projection(self):
        projection = Projection()
        projection.get_projection(self.files[0], self.opt)
        return projection


def main():
    #
    # Read the command line options
    #
    opt, args = get_options()
    #
    # Set debug level
    #
    if opt.quiet:
        print "Extract and join is running in quiet mode"
        log.setLevel("INFO")
    else:
        log.setLevel("DEBUG")

    extract_and_join = ExtractAndJoin(opt, args)
    extract_and_join.run()
    print "Extract and join finished successfully"