import logging
from wrfncxnj import ExtractAndJoin

log = logging.getLogger(__name__)
log.setLevel("DEBUG")


class Options(object):
    pass


def test_extract_and_join():
    opt = Options()
    opt.requested_variables = "T2"
    opt.OFILE = "/tmp/test_cf.nc"
    opt.OUTPUT_PATTERN = None
    opt.TEMPDIR = None
    opt.globfiles = None
    opt.filelist = None
    opt.geofile = None
    opt.fullfile = None
    opt.refdate = None
    opt.splitvars = None
    opt.selected_plevs = None
    opt.vtable = '/home/users/garciam/git/WRFtoolbox_github/WRFToolbox/wrfncxnj/wrfncxnj.table'
    opt.time_units = "hours since 1950-01-01 00:00:00"
    opt.oformat = "NETCDF4_CLASSIC"
    opt.zaxis = None
    opt.paxis = None
    opt.saxis = None
    opt.maxis = None
    opt.tbounds = None
    opt.attributes = None
    opt.discard = None
    opt.singlerec = None
    opt.quiet = None
    opt.ftimes = None
    args = ["/home/users/garciam/pruebas/wrf2cf/wrfout.nc", ]
    extract_and_join = ExtractAndJoin(opt, args)
    extract_and_join.run()


if __name__ == "__main__":
    test_extract_and_join()
