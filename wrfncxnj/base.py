#
# Basic classes and functions use to build wrfncxnj workflow
#
import netCDF4 as ncdf
import numpy as np
import logging
from datetime import datetime
import sys, time, csv, os, shutil, re

log = logging.getLogger(__name__)


class Constants:
    Rd = 287.04
    Rv = 461.5
    RdRv = Rd / Rv
    cp = 7.*Rd/2.
    epsilon_gamma = 0.62197
    es_base_bolton = 0.6112
    es_Abolton = 17.67
    es_Bbolton = 243.5
    es_base_tetens = 6.1078
    es_Atetens_vapor = 7.5
    es_Btetens_vapor = 237.3
    es_Atetens_ice = 9.5
    es_Btetens_ice = 265.5
    g = 9.81
    p1000mb = 100000.
    rcp = Rd/cp
    tkelvin = 273.15


def get_oncvar(
        varobj,
        incvar,
        onc,
        options,
        out_is_2D_but_in_3D = False,
        screenvar_at_2m = False,
        screenvar_at_10m = False,
        screenvar_at_100m = False
):
    oncvar = Oncvar(
        varobj,
        incvar,
        onc,
        options,
        out_is_2D_but_in_3D=out_is_2D_but_in_3D,
        screenvar_at_2m=screenvar_at_2m,
        screenvar_at_10m=screenvar_at_10m,
        screenvar_at_100m=screenvar_at_100m
    )
    return oncvar


def compute_mslp(p, pb, ph, phb, t, qvapor):
    """
    Pure python code by J. Fernandez to extrapolate surface pressure to sea
    level. Strategy borrowed from from_wrf_to_grads.f90 code.
    """
    # Some required physical constants:
    Rd=287.04
    g=9.81
    gamma=0.0065
    # Specific constants for assumptions made in this routine:
    TC=273.16+17.5
    pconst = 10000
    cp = 7.*Rd/2.
    rcp = Rd/cp
    p1000mb = 100000.
    # Transpose and get full variables out of perturbations and potential T
    p = np.transpose(p + pb)
    ph = np.transpose((ph + phb) / 9.81)
    qvapor = np.transpose(qvapor)
    t = np.transpose(t)
    t = (t+300.)*(p/p1000mb)**rcp
    # populate the geopotential_height at mid_levels array with
    # averages between layers below and above
    nz = ph.shape[2]
    z = (ph[:,:,:nz-1] + ph[:,:,1:nz]) / 2.0
    # Find least zeta level that is pconst Pa above the surface. We later use
    # this level to extrapolate a surface pressure and temperature, which is
    # supposed to reduce the effect of the diurnal heating cycle in the pressure
    # field.
    dp = p-(p[:, :, 0, :][:, :, np.newaxis, :]-pconst)
    level = np.add.reduce(dp > 0, axis=2)
    # Get temperature pconst Pa above surface.	Use this to extrapolate
    # the temperature at the surface and down to sea level.
    indic = np.ones(p.shape)*level[:, :, np.newaxis, :]
    loidx = np.arange(nz-1)[np.newaxis, np.newaxis, :, np.newaxis] == indic
    hiidx = np.arange(nz-1)[np.newaxis, np.newaxis, :, np.newaxis] == (indic+1)

    plo = np.add.reduce(p*loidx, 2)
    phi = np.add.reduce(p*hiidx, 2)
    qlo = np.add.reduce(qvapor*loidx, 2)
    qhi = np.add.reduce(qvapor*hiidx, 2)
    tlo = np.add.reduce(t*loidx, 2)*(1. + 0.608 * qlo)
    thi = np.add.reduce(t*hiidx, 2)*(1. + 0.608 * qhi)
    zlo = np.add.reduce(z*loidx, 2)
    zhi = np.add.reduce(z*hiidx, 2)

    p_at_pconst = p[..., 0, :] - pconst
    t_at_pconst = thi-(thi-tlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)
    z_at_pconst = zhi-(zhi-zlo)*np.log(p_at_pconst/phi)*np.log(plo/phi)

    t_surf = t_at_pconst*(p[..., 0, :]/p_at_pconst)**(gamma*Rd/g)
    t_sea_level = t_at_pconst+gamma*z_at_pconst
    # If we follow a traditional computation, there is a correction to the sea
    # level temperature if both the surface and sea level temperatures are
    # *too* hot.
    t_sea_level = np.where((t_sea_level > TC)*(t_surf <= TC),
                           TC, TC - 0.005*(t_surf-TC)**2)
    z_half_lowest = z[:, :, 0, :]
    exponent = (2.*g*z_half_lowest)/(Rd*(t_sea_level + t_surf))
    sea_level_pressure = p[:, :, 0, :] * np.exp(exponent)
    return np.transpose(sea_level_pressure)


def charr2str(carr):
    """
    Convert char (byte) array to string
    """
    return "".join([x.decode('utf8') for x in carr])


def strptime(istr):
    """
    Improved version of strftime that understand any UNIDATA compliant date
    without requiring specifying the format.
    """
    date = ncdf.num2date(0, units="days since %s" % istr)
    return date


def str2offset(istr, time_units):
    if istr == "0000-00-00_00:00:00":
        rval = 0
    else:
        thisdate = strptime(istr)
        rval = ncdf.date2num(thisdate, units=time_units)
    return rval


def discard_suspect_files(filelist, criteria='uncommon_size'):
    total_items = len(filelist)
    file_sizes = map(os.path.getsize, filelist)
    sizes = {}
    for size in file_sizes:
        try:
            sizes[size] += 1
        except:
            sizes[size] = 1
    lsizes = [(b/float(total_items),a) for a,b in sizes.items()]
    lsizes.sort()
    discard_sizes = []
    for item in lsizes:
        if item[0] < 0.2:
            discard_sizes.append(item[1])
    i_file = 0
    rval = []
    for file in filelist:
        print(os.path.basename(file), file_sizes[i_file], end=' ')
        if file_sizes[i_file] in discard_sizes:
            print(" X")
        else:
            rval.append(filelist[i_file])
            print()
        i_file += 1
    return rval


def create_bare_curvilinear_CF_from_wrfnc(
        filename,
        wrfncfile,
        proj,
        ncdf_format,
        time_units,
        options,
        createz=None,
        createp=None,
        createsoil=None,
        createm=None,
        single_level=None
):
    log.debug("Creating %s netCDF file" % filename)
    opt = options
    inc = ncdf.Dataset(wrfncfile, 'r')
    onc = ncdf.Dataset(filename, "w", format=ncdf_format)
    onc.history = "Created by {} on {}".format(sys.argv[0], time.ctime(time.time()))
    onc.sync()
    proj.set_projection(onc)
    if createz:
        if single_level:
            onc.createDimension("z", 1)
        else:
            onc.createDimension("z", len(inc.dimensions["bottom_top"]))
        oncz = onc.createVariable("z", np.float64, ("z",))
        oncz.axis = "Z"
        oncz.long_name = "sigma at layer midpoints"
        oncz.positive = "down"
        oncz.standard_name = "atmosphere_sigma_coordinate"
        oncz.formula_terms = "sigma: z ps: ps ptop: PTOP"
        if "ZNU" in inc.variables and not single_level:
            if len(inc.variables["ZNU"].shape) == 1:
                oncz[:] = inc.variables["ZNU"][:]
            else:
                oncz[:] = inc.variables["ZNU"][0, :]
        if single_level:
            oncz[:] = single_level

        if "P_TOP" in inc.variables:
            oncptop = onc.createVariable("PTOP", np.float32, ())
            oncptop.long_name = "Pressure at the top of the atmosphere"
            oncptop.units = "Pa"
            oncptop[0] = inc.variables["P_TOP"][0]
        onc.sync()
    if createp:
        if single_level:
            onc.createDimension("plev", 1)
        else:
            onc.createDimension("plev",
                                len(inc.dimensions["num_metgrid_levels"]))
        oncp = onc.createVariable("plev", np.float64, ("plev",))
        oncp.axis = "Z"
        oncp.units = "Pa"
        oncp.long_name = "Pressure levels"
        oncp.positive = "down"
        oncp.standard_name = "air_pressure"
        if "PLEV" in inc.variables and not single_level:
            oncp[:] = inc.variables["PLEV"][:]
            #
        # Looks also for the "pressure" variable in hPa so it also works with
        # the original p_interp of NCAR.
        #
        elif "pressure" in inc.variables and not single_level:
            oncp[:] = inc.variables["pressure"][:]*100
        if single_level:
            oncp[:] = single_level
        onc.sync()
    if createsoil:
        if opt.fullfile:
            thisinc = ncdf.Dataset(opt.fullfile,'r')
        else:
            thisinc = inc
        if single_level:
            onc.createDimension("slev", 1)
        else:
            onc.createDimension("slev",
                                len(thisinc.dimensions["soil_layers_stag"]))
        oncsl = onc.createVariable("slev", np.float64, ("slev",))
        oncsl.axis = "Z"
        oncsl.long_name = "Soil level"
        oncsl.units = "m"
        oncsl.positive = "down"
        oncsl.standard_name = "depth_below_surface"
        if "ZS" in thisinc.variables and not single_level:
            oncsl[:] = thisinc.variables["ZS"][0]
        else:
            oncsl[:] = single_level
        onc.sync()
    if createm:
        if single_level:
            onc.createDimension("height", 1)
        else:
            onc.createDimension("height",
                                len(inc.dimensions["num_height_levels"]))
        oncm = onc.createVariable("height", np.float64, ("height",))
        oncm.axis = "Z"
        oncm.units = "m"
        oncm.long_name = "height above the ground"
        oncm.positive = "up"
        oncm.standard_name = "height"
        if "HEIGHT" in inc.variables and not single_level:
            oncm[:] = inc.variables["HEIGHT"][:]
        else:
            oncm[:] = single_level
        onc.sync()
    #
    # Lat-lons (from geo_em or full file if provided)
    #
    if opt.geofile:
        incgeo = ncdf.Dataset(opt.geofile, 'r')
        lats = incgeo.variables["XLAT_M"][0]
        lons = incgeo.variables["XLONG_M"][0]
        incgeo.close()
    elif opt.fullfile:
        incgeo = ncdf.Dataset(opt.fullfile, 'r')
        lats = incgeo.variables["XLAT"][0]
        lons = incgeo.variables["XLONG"][0]
        incgeo.close()
    else:
        lats = inc.variables["XLAT"][0]
        lons = inc.variables["XLONG"][0]
    onclat = onc.createVariable("lat", np.float64, (proj.yname, proj.xname))
    onclat.long_name = "latitude"
    onclat.standard_name = "latitude"
    onclat.units = "degrees_north"

    onclat[:len(lats)] = lats
    onclon = onc.createVariable("lon",np.float64, (proj.yname, proj.xname))
    onclon.long_name = "longitude"
    onclon.standard_name = "longitude"
    onclon.units = "degrees_east"
    onclon[:len(lons)] = lons
    #
    # Get the initial date and create a new time variable
    #
    onc.createDimension("time", None)
    onctime = onc.createVariable("time", np.float64, ("time",))
    onctime.long_name = "time"
    onctime.standard_name = "time"
    onctime.units = time_units
    onctime.calendar = "standard"
    if opt.tbounds:
        onctime.bounds = "time_bnds"
        onc.createDimension("nv", 2)
        onc.createVariable("time_bnds", np.float64, ("time", "nv"))
    inc.close()
    onc.sync()

    return onc


def add_height_coordinate(onc, coorname, val):
    if coorname not in onc.variables:
        onc.createDimension(coorname,1)
        #hvar = onc.createVariable(coorname, 'f', (coorname,))
        hvar = onc.createVariable(coorname, 'float64', (coorname,))
        hvar.long_name = "height above the ground"
        hvar.standard_name = "height"
        hvar._CoordinateAxisType = "Height"
        hvar.units = "m"
        hvar[0] = np.array(val, 'f')


def add_depth_coordinate(onc, coorname, val):
    # TODO: Hay que aniadir las boundaries...
    if coorname not in onc.variables:
        onc.createDimension(coorname, 1)
        hvar = onc.createVariable(coorname, 'f', (coorname,))
        hvar.long_name = "depth below the surface"
        hvar.standard_name = "depth"
        hvar.units = "m"
        hvar[0] = np.array(val, 'f')


projection_map = {
    1: "Lambert_Conformal",
    2: "Polar_Stereographic",
    3: "Mercator",
    6: "Rotated_Pole"
}


class Projection:
    def get_projection(self, wrfncfile, options):
        wrfncobj = ncdf.Dataset(wrfncfile, "r")
        self.projection = projection_map[wrfncobj.MAP_PROJ]
        if self.projection == "Lambert_Conformal":
            self.grid_mapping_name = "lambert_conformal_conic"
            self.cone_type = "secant"
            self.northern_parallel = "%4.1f" % wrfncobj.TRUELAT2
            self.southern_parallel = "%4.1f" % wrfncobj.TRUELAT1
            self.longitude_of_central_meridian = "%4.1f" % wrfncobj.CEN_LON
            self.latitude_of_projection_origin = "%4.1f" % wrfncobj.CEN_LAT
            self.proj_attr = [
                "grid_mapping_name",
                "cone_type",
                "northern_parallel",
                "southern_parallel",
                "longitude_of_central_meridian",
                "latitude_of_projection_origin"
            ]
            # CF attributes for x and y coordinates.
            self.xname = "x"
            self.xlen = len(wrfncobj.dimensions["west_east"])
            self.xvalues = (np.arange(1, self.xlen + 1) - self.xlen/2)*wrfncobj.DX/1000
            self.xattr = {
                "standard_name" : "projection_x_coordinate",
                "long_name" : "x coordinate of projection",
                "axis" : "X"	,
                "units" : "km"
            }
            self.yname = "y"
            self.ylen = len(wrfncobj.dimensions["south_north"])
            self.yvalues = (np.arange(1, self.ylen + 1) - self.ylen/2)*wrfncobj.DY/1000
            self.yattr = {
                "standard_name" : "projection_y_coordinate",
                "long_name" : "y coordinate of projection",
                "axis" : "Y",
                "units" : "km"
            }
        elif self.projection == "Rotated_Pole":
            self.grid_mapping_name = "rotated_latitude_longitude"
            self.grid_north_pole_latitude = "%4.1f" % wrfncobj.POLE_LAT
            self.grid_north_pole_longitude = "%4.1f" % wrfncobj.POLE_LON
            self.proj_attr = [
                "grid_mapping_name",
                "grid_north_pole_latitude",
                "grid_north_pole_longitude"
            ]
            # CF attributes for x and y coordinates.
            self.xname = "rlon"
            self.xlen = len(wrfncobj.dimensions["west_east"])
            self.xattr = {
                "standard_name" : "grid_longitude",
                "long_name" : "longitude in rotated pole grid",
                "axis" : "X",
                "units" : "degrees"
            }
            self.yname = "rlat"
            self.ylen = len(wrfncobj.dimensions["south_north"])
            self.yattr = {
                "standard_name" : "grid_latitude",
                "long_name" : "latitude in rotated pole grid",
                "axis" : "Y",
                "units" : "degrees"
            }
            if options.geofile:
                incgeo = ncdf.Dataset(options.geofile,'r')
                self.xvalues = incgeo.variables["CLONG"][0, 0, :]
                self.yvalues = incgeo.variables["CLAT"][0, :, 0]
                incgeo.close()
            else:
                self.xvalues = wrfncobj.variables["CLONG"][0, 0, :]
                self.yvalues = wrfncobj.variables["CLAT"][0, :, 0]
        elif self.projection == "Mercator":
            self.grid_mapping_name = "mercator"
            self.longitude_of_projection_origin = "%4.1f" % wrfncobj.CEN_LON
            self.standard_parallel = "%4.1f" % wrfncobj.TRUELAT1
            self.proj_attr = [
                "grid_mapping_name",
                "longitude_of_projection_origin",
                "standard_parallel"
            ]
            # CF attributes for x and y coordinates.
            self.xname = "x"
            self.xlen = len(wrfncobj.dimensions["west_east"])
            self.xvalues = (np.arange(1, self.xlen + 1) - self.xlen/2)*wrfncobj.DX/1000
            self.xattr = {
                "standard_name" : "projection_x_coordinate",
                "long_name" : "x coordinate of projection",
                "axis" : "X"	,
                "units" : "km"
            }
            self.yname = "y"
            self.ylen = len(wrfncobj.dimensions["south_north"])
            self.yvalues = (np.arange(1, self.ylen + 1) - self.ylen/2)*wrfncobj.DY/1000
            self.yattr = {
                "standard_name" : "projection_y_coordinate",
                "long_name" : "y coordinate of projection",
                "axis" : "Y",
                "units" : "km"
            }
        else:
            raise NotImplementedError("Error: Projection %s is currently not "
                                      "supported by XnJ" % self.projection)
        wrfncobj.close()
        #
        # Dictionary for Oncvar.__init__
        #
        self.dimension_map = {
            "south_north": self.yname,
            "west_east": self.xname,
            "bottom_top": "z",
            "soil_layers_stag": "slev",
            "num_metgrid_levels": "plev",
            "bottom_top_stag": "z", # the variables should be de-staggered before copying them to the output file
            "num_height_levels" : "height"
        }

    def set_projection(self, onc):
        oncproj = onc.createVariable(self.projection, "c", ())
        for attr in self.proj_attr:
            setattr(oncproj, attr, getattr(self, attr))
        onc.createDimension(self.xname, self.xlen)
        onc.createDimension(self.yname, self.ylen)
        oncx = onc.createVariable(self.xname, np.float64, (self.xname,))
        oncy = onc.createVariable(self.yname, np.float64, (self.yname,))
        for attr, attrval in self.xattr.items():
            setattr(oncx, attr, attrval)
        for attr, attrval in self.yattr.items():
            setattr(oncy, attr, attrval)
        oncx[:] = self.xvalues
        oncy[:] = self.yvalues
        onc.sync()


class Variable:
    def __repr__(self):
        return """\
abbr: %(standard_abbr)s
standard_name: %(standard_name)s
units: %(units)s
""" % self.__dict__

    def get_levels(self, wrfncfile, options):
        opt = options
        inc = ncdf.Dataset(wrfncfile, 'r')
        #
        # Soil levels
        #
        if self.saxis:
            if "ZS" in inc.variables:
                levels = inc.variables["ZS"][0, :]
            elif opt.fullfile:
                incfull = ncdf.Dataset(opt.fullfile,'r')
                levels = incfull.variables["ZS"][0, :]
                incfull.close()
        #
        # Pressure levels
        #
        elif self.paxis:
            if "PLEV" in inc.variables:
                levels = inc.variables["PLEV"][:]
            #
            # Looks also for the "pressure" variable in hPa so it also works
            # with the original p_interp of NCAR.
            #
            elif "pressure" in inc.variables:
                levels = inc.variables["pressure"][:]*100
            else:
                raise RuntimeError(
                    "Error, no pressure levels found in the input files")
        #
        # Eta levels
        #
        elif self.zaxis:
            if "ZNU" in inc.variables:
                if len(inc.variables["ZNU"].shape) == 1:
                    levels = inc.variables["ZNU"][:]
                else:
                    levels = inc.variables["ZNU"][0]
            elif opt.fullfile:
                incfull = ncdf.Dataset(opt.fullfile, 'r')
                if len(incfull.variables["ZNU"].shape) == 1:
                    levels = incfull.variables["ZNU"][:]
                else:
                    levels = incfull.variables["ZNU"][0]
                incfull.close()
            else:
                raise RuntimeError(
                    "Error, no eta levels found in the input files")
        #
        # Height levels in meters
        #
        elif self.maxis:
            if "HEIGHT" in inc.variables:
                levels = inc.variables["HEIGHT"][:]
        else:
            raise RuntimeError(
                "Error, the var %s has not any vertical axis "
                "defined in the wrfncxnj_table." % self.abbr)
        inc.close()
        self.levels = levels
        self.level_map = {}
        for i in range(0, len(levels)):
            self.level_map[levels[i]] = i
        if opt.selected_plevs and self.paxis:
            for lev in levels:
                if lev not in np.array(opt.selected_plevs.split(","),
                                       dtype="float")*100:
                    del self.level_map[lev]
        if opt.selected_slevs and self.saxis:
            for lev in levels:
                slev = "%.2f" % lev
                if slev not in opt.selected_slevs.split(","):
                    del self.level_map[lev]


class WrfNcFiles:
    def __init__(self):
        self.current = None
        self.nxt = None
        self.prv = None
        self.geo = None
        self.full = None

    def loadFiles(self, filelist, previous_file=None, next_file=None):
        self.filelist = filelist
        self.nfiles = len(filelist)
        self.next_file = next_file
        self.ifile = -1
        if previous_file:
            self.current = ncdf.Dataset(previous_file, "r")
        else:
            self.current = None
        self.nxt = ncdf.Dataset(self.filelist[0])

    def assignNext(self):
        if self.ifile == self.nfiles-1:
            if self.next_file:
                self.nxt = ncdf.Dataset(self.next_file, "r")
            else:
                self.nxt = None
        else:
            self.nxt = ncdf.Dataset(self.filelist[self.ifile+1], "r")

    def cycle(self):
        self.ifile += 1
        if self.prv:
            self.prv.close()
        self.prv = self.current
        self.current = self.nxt
        self.assignNext()

    def __iter__(self):
        return self

    def __next__(self):
        if self.ifile >= self.nfiles-1:
            raise StopIteration
        else:
            self.cycle()
            return self

    def rewind(self):
        self.loadFiles(self.filelist)

    def close(self):
        if self.current:
            self.current.close()
        if self.nxt:
            self.nxt.close()
        if self.prv:
            self.prv.close()

    def closeCommon(self):
        if self.geo:
            self.geo.close()
        if self.full:
            self.full.close()

    def __str__(self):
        return """
            ifile: %(ifile)i
            prev: %(prv)s
            curr: %(current)s
            next: %(nxt)s
        """ % self.__dict__


class WrfNcTime:
    def __init__(self, initialdate, time_units):
        self.is_singlerec = False
        self.nrec = 0
        self.initialdate = initialdate
        self.iini = 0
        self.iend = None
        self.outstep_s = 0
        self.time_units = time_units

    def checkStep(self, wrfnc):
        # TODO: Check here that the next file does not have a time gap
        incTimes = wrfnc.current.variables["Times"]
        t0 = strptime(charr2str(incTimes[0]))
        if len(incTimes) > 1:
            t1 = strptime(charr2str(incTimes[1]))
            delta = t1-t0
        elif wrfnc.nxt:
            t1 = strptime(
                charr2str(wrfnc.nxt.variables["Times"][0]))
            delta = t1-t0
        elif wrfnc.prv:
            tp = strptime(
                charr2str(wrfnc.prv.variables["Times"][:][-1]))
            delta = t0-tp
        else:
            delta = t0 - t0
        self.outstep_s = 86400*delta.days + delta.seconds

    def getTimes(self, wrfnc, is_singlerec):
        incTimes = wrfnc.current.variables["Times"]
        self.nrec = len(incTimes)
        if is_singlerec:
            self.nrec = 1
        self.iend = self.iini + self.nrec
        times = map(charr2str, incTimes[:self.nrec])
        return [str2offset(x, self.time_units) for x in times]

    def cycle(self):
        self.iini += self.nrec

    def __str__(self):
        return """
                iini = %(iini)d
                iend = %(iend)d
                nrec = %(nrec)d
                outstep_s = %(outstep_s)d
        """ % self.__dict__


def read_variable_table(vars, vtable, proj, wrfncfile, options):
    rval = {}
    csv_reader = csv.reader(open(vtable), delimiter=";",
                            skipinitialspace=True)
    for line in csv_reader:
        if line[0].strip() == "":
            continue
        if line[0][0] == "#":
            continue
        if len(line) > 5 and line[5] == "do_nothing":
            continue
        varwrf = line[0].strip()
        if varwrf in vars:
            v = Variable()
            v.varname = varwrf
            v.standard_abbr = line[1].strip()
            v.long_name = line[2].strip()
            v.standard_name = line[3].strip()
            v.units = line[4].strip()
            try:
                v.transform = line[5].strip()
            except:
                v.transform = ""
            try:
                v.paxis = line[6].strip() == "p"
                v.zaxis = line[6].strip() == "z"
                v.saxis = line[6].strip() == "s"
                v.maxis = line[6].strip() == "m"
            except:
                v.paxis = v.zaxis = v.saxis = v.maxis = False
            v.is_3D = v.paxis or v.zaxis or v.saxis or v.maxis
            v.projection = proj.projection
            v.dimension_map = proj.dimension_map
            if v.is_3D:
                v.get_levels(wrfncfile, options)
            rval[varwrf] = v
    return rval


class Oncvar:
    def __init__(
            self,
            varobj,
            incvar,
            onc,
            options,
            out_is_2D_but_in_3D=False,
            screenvar_at_2m=False,
            screenvar_at_10m=False,
            screenvar_at_100m=False
    ):
        self.screenvar_at_2m = screenvar_at_2m
        self.screenvar_at_10m = screenvar_at_10m
        self.screenvar_at_100m = screenvar_at_100m
        self.incvar = incvar
        self.options = options
        if out_is_2D_but_in_3D:
            cut_from = 2
        else:
            cut_from = 1
        variable_dims_except_time = tuple(map(lambda x: varobj.dimension_map[x],
                                              incvar.dimensions[cut_from:]))
        self.dims = ("time",) + variable_dims_except_time
        if screenvar_at_2m:
            self.dims = self.dims[:1] + ("height",) + self.dims[1:]
        if screenvar_at_10m:
            self.dims = self.dims[:1] + ("heightv",) + self.dims[1:]
        self.onc = onc
        self.dtype = np.float32

    def set_ncvars(self, varobj):
        opt = self.options
        self.varobj = varobj
        #
        # If varobj is a 3D variable and we want to split it, then self.ncvars
        # is a dict, else it is a single ncdf.Variable object.
        #
        if varobj.is_3D and opt.splitlevs:
            self.ncvars = {}
            if varobj.standard_abbr in self.onc[varobj.level_map.keys()[0]].variables:
                for level in np.sort(self.onc.keys()):
                    self.ncvars[level] = self.onc[level].variables[varobj.standard_abbr]
            else:
                for level in np.sort(self.onc.keys()):
                    if opt.oformat == "NETCDF4_CLASSIC":
                        self.ncvars[level] = self.onc[level].createVariable(
                                varobj.standard_abbr,
                                np.float32,
                                self.dims,
                                zlib=True,
                                complevel=4,
                                shuffle=True
                            )
                    else:
                        self.ncvars[level] = self.onc[level].createVariable(
                            varobj.standard_abbr,
                            np.float32,
                            self.dims
                        )
                    set_ncvar_attributes(self.ncvars[level], varobj)
        else:
            if varobj.standard_abbr in self.onc.variables:
                self.ncvars = self.onc.variables[varobj.standard_abbr]
            else:
                if self.screenvar_at_2m:
                    add_height_coordinate(self.onc, "height", 2)
                if self.screenvar_at_10m:
                    add_height_coordinate(self.onc, "heightv", 10)
                    if self.screenvar_at_100m:
                            add_height_coordinate(self.onc, "height", 100)
                if opt.oformat == "NETCDF4_CLASSIC":
                    self.ncvars = self.onc.createVariable(
                        varobj.standard_abbr,
                        np.float32,
                        self.dims,
                        zlib=True,
                        complevel=4,
                        shuffle=True
                    )
                else:
                    self.ncvars = self.onc.createVariable(
                        varobj.standard_abbr,
                        np.float32,
                        self.dims
                    )
                set_ncvar_attributes(self.ncvars, varobj)
    #
    # If splitting levels, second dimension of the array is assumed to be
    # the vertical coordinate.
    #
    def assign_values(self, array, ini, end, times):
        opt = self.options
        if self.varobj.is_3D:
            if not len(self.varobj.levels) == array.shape[1]:
                raise RuntimeError(
                    "Error: Input variable %s has a wrong number of "
                    "levels." % self.incvar
                )
        if self.varobj.is_3D and opt.splitlevs:
            for level in np.sort(self.varobj.level_map.keys()):
                i = self.varobj.level_map[level] # get the level index
                log.debug("i = {} Filling level {}".format(i, level))
                self.ncvars[level][ini:end, :] = array[:, i, :, :]
                self.onc[level].variables["time"][ini:end] = times
        else:
            self.ncvars[ini:end, :] = array[:]
            self.onc.variables["time"][ini:end] = times

    def sync(self):
        if self.varobj.is_3D and self.options.splitlevs:
            for nc in self.onc.values():
                nc.sync()
        else:
            self.onc.sync()

    def close(self):
        if self.varobj.is_3D and self.options.splitlevs:
            for nc in self.onc.values():
                nc.close()
        else:
            self.onc.close()

    def __repr__(self):
        return """
dims: %(dims)s
onc: %(onc)s
ncvars: %(ncvars)s
""" % self.__dict__


def map_levels(varobj):
    levelmap = {}
    for i, level in enumerate(varobj.levels[::-1]):
        levelmap[i] = level
    return levelmap


def set_ncvar_attributes(ncvar, varobj):
    ncvar.long_name = varobj.long_name
    ncvar.standard_name = varobj.standard_name
    ncvar.units = varobj.units
    ncvar.coordinates="lat lon"
    ncvar.grid_mapping = varobj.projection


def get_first_and_last_times(ifiles, format):
    firstnc = ncdf.Dataset(ifiles[0], "r")
    lastnc = ncdf.Dataset(ifiles[-1], "r")
    ftimes = firstnc.variables["Times"][:]
    lasttimes = lastnc.variables["Times"][:]
    firsttime = strptime(charr2str(ftimes[0]))
    lasttime = strptime(charr2str(lasttimes[-1]))
    firstnc.close()
    lastnc.close()
    return firsttime, lasttime


def replace_output_pattern(varobj, pattern, firstdate, lastdate, options,
                           level=None):
    """
    Parse file name patterns replacing [varcf] and [varwrf]
    """
    pout = pattern.replace('[varcf]', varobj.standard_abbr)
    pout = pout.replace('[varwrf]', varobj.varname)
    if options.ftimes:
        pout = pout.replace('[firsttime]', options.ftimes.split(",")[0])
        pout = pout.replace('[lasttime]', options.ftimes.split(",")[1])
    elif options.singlerec:
        pass
    else:
        pout = pout.replace('[firsttime]', firstdate.strftime("%Y%m%d%H"))
        pout = pout.replace('[lasttime]', lastdate.strftime("%Y%m%d%H"))
    if level:
        if varobj.paxis:
            # Assuming that the level is in Pa and that we want it in the file
            # name in hPa, to save space.
            pout = pout.replace('[level]', "%g" % (float(level)/100.))
        else:
            pout = pout.replace('[level]', "%g" % level)
    else:
        pout = pout.replace('[level]', "")
    return pout
#
# Set global attributes.
#
def set_global_attributes(ncobj, options):
    gattr_file = csv.reader(open(options.attributes), delimiter=" ",
                            skipinitialspace=True)
    for line in gattr_file:
        setattr(ncobj, line[0], line[1])


def get_onc_lon_lat_time(oncs, options):
    if options.splitvars:
        try:
            outvars = oncs[oncs.keys()[0]].variables
            outnc = oncs[oncs.keys()[0]]
        except:
            outncdict = oncs[oncs.keys()[0]]
            outnc = outncdict[outncdict.keys()[0]]
    else:
        outnc = oncs
    lon = outnc.variables["lon"]
    lat = outnc.variables["lat"]
    time = outnc.variables["time"]
    return lon, lat, time


class OfileHandler:
    """
    Handle output files and paths.
    """
    def __init__(self, options):
        self.options = options
        opt = options
        if opt.OFILE:
            self.ofile = opt.OFILE
        if opt.OUTPUT_PATTERN:
            self.ofile = opt.OUTPUT_PATTERN
        if not os.path.dirname(self.ofile):
            self.ofile = os.getcwd() + "/" + self.ofile
            self.odir = os.path.dirname(self.ofile)
        if opt.TEMPDIR:
            self.tempdir = opt.TEMPDIR
            log.debug("Writing %s files using %s as temporary directory" %
                      (self.ofile, self.tempdir))

    def set_oncnames(self, oncnames):
        self.oncnames = oncnames

    def set_tempdir(self):
        opt = self.options
        if opt.TEMPDIR:
            if os.path.exists(self.tempdir):
                print("Warning: %s already exists, overwriting" % self.tempdir)
                shutil.rmtree(self.tempdir)
            os.makedirs(self.tempdir)
            self.tempfile = self.ofile.replace(os.path.dirname(self.ofile),
                                               self.tempdir)
        else:
            self.tempfile = self.ofile

    def filter_times(self, itime, ftime):
        for file in self.oncnames:
            ofile_filtered = file.replace(".nc", "_filtered.nc")
            netcdf_seldate(itime, ftime, file, ofile_filtered, self.options)
            os.remove(file)
            os.rename(ofile_filtered, file)

    def move_output(self):
        for file in self.oncnames:
            shutil.move(file, os.path.dirname(os.path.abspath(self.ofile)))
        shutil.rmtree(self.tempdir)


def ofile_exists(ofile, ofh, options):
    if options.TEMPDIR:
        ofile_finaldir = os.path.dirname(ofh.ofile) + "/"  + os.path.basename(ofile)
    else:
        ofile_finaldir = ofile
    return os.path.exists(ofile_finaldir)


def create_oncs(vars, ofh, ifiles, wnt, proj, options):
    oncs = {}
    oncnames = []
    if options.splitvars:
        if options.singlerec:
            firstt, lastt = None, None
        else:
            firstt, lastt = get_first_and_last_times(ifiles, format)
        for var in vars.values():
            log.debug(var.varname)
            if options.splitlevs and var.is_3D:
                log.debug(np.sort(var.level_map.keys()))
                # Define a dict of oncs for the 3D variables.
                oncs[var.varname] = {}
                for level in np.sort(var.level_map.keys()):
                    ofile = replace_output_pattern(var, ofh.tempfile, firstt,
                                                   lastt, options, level)
                    # Check if the file does already exists in the output folder.
                    if ofile_exists(ofile, ofh, options):
                        log.warning("File %s already exists, skipping "
                                    "variable %s, "
                                    "level %s" % (ofile, var.varname, level)
                                    )
                        del var.level_map[level]
                        # Check if the level_map dic is now empty.
                        if not var.level_map.keys():
                            del vars[var.varname]
                            del oncs[var.varname]
                            if not vars.keys():
                                log.info("All the files have been already "
                                         "processed, exiting")
                                sys.exit(0)
                        continue
                    if os.path.exists(ofile):
                        errmsg = ("Error: File %s has been already created in "
                                  "the temporal folder. Change --output-pattern"
                                  " or variable list so output file names do"
                                  " not overlap." % (ofile)
                                  )
                        raise OSError(errmsg)
                    cf_netcdf = create_bare_curvilinear_CF_from_wrfnc(
                        ofile,
                        ifiles[0],
                        proj,
                        options.oformat, # TODO: remove and read from options
                        options.time_units,
                        options,
                        var.zaxis,
                        var.paxis,
                        var.saxis,
                        var.maxis,
                        level
                    )
                    oncs[var.varname][level] = cf_netcdf
                    if options.attributes:
                        set_global_attributes(oncs[var.varname][level],options)
                    oncnames.append(ofile)
            else:
                ofile = replace_output_pattern(var, ofh.tempfile, firstt,
                                               lastt, options)
                # Check if the file does already exists.
                if ofile_exists(ofile, ofh, options):
                    log.warning("Warning: File %s already exists, skipping "
                                "variable %s" % (ofile, var.varname))
                    del vars[var.varname]
                    if not vars.keys():
                        print("All the files have been already processed, exiting")
                        sys.exit(0)
                    continue
                oncs[var.varname] = create_bare_curvilinear_CF_from_wrfnc(
                    ofile,
                    ifiles[0],
                    proj,
                    options.oformat,
                    options.time_units,
                    options,
                    var.zaxis,
                    var.paxis,
                    var.saxis,
                    var.maxis
                )
                if options.attributes:
                    set_global_attributes(oncs[var.varname], options)
                oncnames.append(ofile)
    else:
        ofile = options.OFILE
        if ofile_exists(ofile, ofh, options):
            log.warning("File %s already exists, refusing to "
                        "overwrite it" % ofile)
            sys.exit(0)
        oncs = create_bare_curvilinear_CF_from_wrfnc(
            ofh.tempfile,
            ifiles[0],
            proj,
            options.oformat,
            options.time_units,
            options,
            options.zaxis,
            options.paxis,
            options.saxis,
            options.maxis
        )
        if options.attributes:
            set_global_attributes(oncs)
        oncnames.append(ofile)
    if len(set(oncnames)) < len(oncnames):
        log.warning("Two variables are being processed with the same output "
                    "file name, check OUTPUT_PATTERN")
    return oncs, oncnames
#
# Function to copy netCDF structures.
#
def copy_netcdf_structure(ifile, ofile, variables, isncobj = False):
    print("Creating %s netCDF file" % (ofile))
    if isncobj:
        inc = ifile
    else:
        inc = ncdf.Dataset(ifile, 'r')
    onc = ncdf.Dataset(ofile, "w", format = inc.file_format)
    if variables == "all":
        variables = inc.variables.keys()
    #
    # Copy global attributes and redefine history
    #
    for ikey, ivalue in inc.__dict__.iteritems():
        onc.setncattr(ikey, ivalue)
    onc.history = "Created by {} on {}".format(sys.argv[0],time.ctime(time.time()))
    onc.sync()
    #
    # Copy dimensions
    #
    for dimname, dimobj in inc.dimensions.iteritems():
        print("Setting dimension {} {}".format(dimname, dimobj))
        if dimobj.isunlimited():
            print("Dimension is unlimited")
            onc.createDimension(dimname, None)
        else:
            onc.createDimension(dimname, len(dimobj))
    onc.sync()
    #
    # Copy variables specified in the argument
    #
    for ivarname, ivarobj in inc.variables.iteritems():
        if ivarname not in variables:
            continue
        else:
            ovarobj = onc.createVariable(ivarname, ivarobj.dtype,
                                         ivarobj.dimensions)
            for attrname, attrvalue in ivarobj.__dict__.iteritems():
                ovarobj.setncattr(attrname, attrvalue)
            if ivarobj.dtype != "S1":
                ovarobj[:] = ivarobj[:]
    #
    # sync and return onc
    #
    onc.sync()
    return onc


def netcdf_seldate(idate, fdate, ifile, ofile, options):
    log.debug("Filtering output files to remove dates out from the interval"
              " %s %s" % (idate, fdate))
    inc = ncdf.Dataset(ifile, "r")
    #
    # Divide the variables in those to directly copy and those that need to
    # be filtered.
    #
    copyvars = inc.variables.keys()
    main_vars = []
    for ivar, ivarobj in inc.variables.iteritems():
        if "time" in ivarobj.dimensions:
            main_vars.append(ivar)
            copyvars.remove(ivar)
    #
    # Filter the main variables to keep the selected times
    #
    t = inc.variables["time"]
    idatetimes = ncdf.num2date(t[:], units=t.units)
    idateobj = datetime.strptime(idate, '%Y%m%d%H')
    fdateobj = datetime.strptime(fdate, '%Y%m%d%H')
    tmask = np.array(
        [ idateobj <= idatetimes[i] <= fdateobj for i in range(0, len(t[:]))])
    main_vars.remove("time")
    #
    # Copy the netCDF structure and record the filtered variables
    #
    onc = copy_netcdf_structure(inc, ofile, variables = copyvars, isncobj=True)
    onct = onc.createVariable("time", t.dtype, t.dimensions)
    for attrname, attrvalue in t.__dict__.iteritems():
        onct.setncattr(attrname, attrvalue)
    onct[:] = ncdf.date2num(idatetimes[tmask], t.units)
    for main_variable in main_vars:
        mainvarobj = inc.variables[main_variable]
        if options.oformat == "NETCDF4_CLASSIC":
            oncvarobj = onc.createVariable(
                main_variable,
                mainvarobj.dtype,
                mainvarobj.dimensions,
                zlib=True, complevel=4,
                shuffle=True
            )
        else:
            oncvarobj = onc.createVariable(
                main_variable,
                mainvarobj.dtype,
                mainvarobj.dimensions
            )
        for attrname, attrvalue in mainvarobj.__dict__.iteritems():
            oncvarobj.setncattr(attrname, attrvalue)
        log.debug("Filtering time axis for {}".format(ivarobj))
        oncvarobj[:] = mainvarobj[tmask, :]
        onc.sync()
