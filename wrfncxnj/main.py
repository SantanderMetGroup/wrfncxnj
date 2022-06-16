import netCDF4 as ncdf
import logging
import os
import re
import string
import sys
from glob import glob
import wrfncxnj.base as base
import wrfncxnj.generic as generic
import wrfncxnj.diagnostics as diagnostics
from wrfncxnj.options import get_options, check_options_consistency

log = logging.getLogger(__name__)


class ExtractAndJoin:
    """
    Class where the main workflow of the program is defined
    """
    def __init__(self, opt, args):
        self.opt = opt
        self.args = args
        self.wrf_files_iterator = base.WrfNcFiles()
        if opt.geofile:
            self.wrf_files_iterator.geo = ncdf.Dataset(opt.geofile, "r")
        if opt.fullfile:
            self.wrf_files_iterator.full = ncdf.Dataset(opt.fullfile, "r")
        check_options_consistency(opt)
        self.ofh = base.OfileHandler(opt)
        self.ofh.set_tempdir()
        self.is_geofile = False
        self.files = self.get_input_files()
        self.wrf_files_iterator.loadFiles(self.files)
        self.requested_variable_names = opt.requested_variables.split(',')
        if not opt.vtable:
            self.opt.vtable = os.path.dirname(__file__) + "/wrfncxnj.table"

        self.projection = self.get_projection()
        self.variables = base.read_variable_table(
            self.requested_variable_names,
            opt.vtable,
            self.projection,
            self.files[0],
            opt
        )
        refdate = opt.time_units.split()[2]
        self.wrfnctime = base.WrfNcTime(base.strptime(refdate), opt.time_units)

    def run(self):
        opt = self.opt
        #
        # Create output netCDF files
        #
        oncs, oncnames = base.create_oncs(
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
            print("Processing var %s" % varname)

        function_name = "compute_%s" % varname
        if hasattr(diagnostics, function_name):
            log.debug("Using the %s function" % function_name)
            process_func = getattr(diagnostics, function_name)
            oncvar, copyval = process_func(varobj, onc, wrfnc, wnt, self.opt)
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
            oncvar = base.get_oncvar(varobj, incvar, onc, self.opt)
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
            files = map(string.strip, open(opt.filelist).readlines())
        else:
            files = self.args
        #
        # Look for geofile
        #
        if not files and self.wrf_files_iterator.geo:
            print("No input files provided.")
            print("Trying to find the variables in the geo_em file provided.")
            files = [opt.geofile, ]
            self.is_geofile = True
        return files

    def get_projection(self):
        projection = base.Projection()
        projection.get_projection(self.files[0], self.opt)
        return projection


class ParseTransform:
    """
    Parses and executes arithmetic expressions in wrfncxnj.table
    """
    def __init__(self, transformstr, options):
        self.transformstr = transformstr
        self.options = options
        words = ""
        lastchar = "X"
        for char in transformstr:
            if char in ["+", "-"] and lastchar == "E":
                # Takes care of constants in exp notation, e.g. 4.8e-3
                words += char
            elif char in ["*", "+", "-", "(", ")", "[", "]", ":", "/"]:
                words += " "
            else:
                words += char
            lastchar = char.upper()
        words = words.split()
        self.variables = []
        self.constants = []
        self.functions = []
        for word in words:
            try:
                self.constants.append(float(word))
            except:
                if word == word.upper():
                    self.variables.append(word)
                else:
                        self.functions.append(word)

    def execute(self, varobj, onc, wnfiles, wntimes, vardic):
        cmdstr = self.transformstr
        if self.variables:
            for var in set(self.variables):
                cmdstr = re.sub(r"\b%s\b" % var, "vardic['%s']" % var, cmdstr)
            log.debug("Executing -> %s" % cmdstr)
            exec("copyval = %s" % cmdstr)
            oncvar = base.get_oncvar(
                varobj,
                wnfiles.current.variables[self.variables[0]],
                onc,
                self.options
            )
        else:
            log.debug("Using function: %s" % self.functions[0])
            process_func = getattr(generic, self.functions[0])
            oncvar, copyval = process_func(varobj, onc, wnfiles, wntimes,
                                           self.options)
        return oncvar, copyval

    def __str__(self):
        return """
            variables: %(variables)s
            functions: %(functions)s
            constants: %(constants)s
        """ % self.__dict__
