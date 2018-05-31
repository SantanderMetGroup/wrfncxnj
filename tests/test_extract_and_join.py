import logging
from wrfncxnj.options import get_options
from wrfncxnj.main import ExtractAndJoin

logging.basicConfig()
log = logging.getLogger()
log.setLevel("DEBUG")


class Options(object):
    pass


def test_ofile():
    opt, noargs = get_options()
    opt.requested_variables = "T2,U10ER,V10ER,RAINF,TEMP"
    opt.paxis = True
    opt.OFILE = "/tmp/test_cf.nc"
    opt.quiet = False
    opt.vtable = '/home/users/garciam/git/WRFtoolbox_github/WRFToolbox/wrfncxnj/wrfncxnj/wrfncxnj.table'
    opt.geofile = "/predictia/Projects/metromadrid/domains/metromad_lakes/geo_em.d03.nc"
    opt.time_units = "hours since 1950-01-01 00:00:00"
    opt.oformat = "NETCDF4_CLASSIC"
    opt.ftimes = "2016070412,2016070418"
    args = ["/home/users/garciam/pruebas/wrf2cf/wrfout.nc", ]
    extract_and_join = ExtractAndJoin(opt, args)
    extract_and_join.run()


def test_output_pattern():
    opt, noargs = get_options()
    opt.requested_variables = "T2,U10ER,V10ER,RAINF,TEMP,GLW"
    opt.splitvars = True
    opt.splitlevs = True
    opt.OUTPUT_PATTERN = "/tmp/[varcf]_[varwrf]_[level]_[firsttime]_[lasttime]_test.nc"
    opt.vtable = '/home/users/garciam/git/WRFtoolbox_github/WRFToolbox/wrfncxnj/wrfncxnj/wrfncxnj.table'
    opt.geofile = "/predictia/Projects/metromadrid/domains/metromad_lakes/geo_em.d03.nc"
    opt.time_units = "hours since 1950-01-01 00:00:00"
    opt.ftimes = "2016070412,2016070418"
    opt.oformat = "NETCDF4_CLASSIC"
    args = ["/home/users/garciam/pruebas/wrf2cf/wrfout.nc", ]
    extract_and_join = ExtractAndJoin(opt, args)
    extract_and_join.run()


if __name__ == "__main__":
    test_ofile()
    #test_output_pattern()
