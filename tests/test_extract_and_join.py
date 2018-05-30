import logging
from wrfncxnj import ExtractAndJoin

log = logging.getLogger(__name__)
log.setLevel("DEBUG")


def test_extract_and_join():
    opt = {
        "requested_variables": ["T2", ],
        "OFILE": "/tmp/test_cf.nc"
    }
    args = ["/home/users/garciam/pruebas/wrf2cf/wrfout.nc", ]
    extract_and_join = ExtractAndJoin(opt, args)
    extract_and_join.run()


if __name__ == "__main__":
    test_extract_and_join()
