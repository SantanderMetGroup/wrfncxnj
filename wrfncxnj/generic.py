# -*- coding: utf-8 -*-
#
#  Generic functions to process WRF variables
#
import numpy as np
from wrfncxnj.base import get_oncvar

#
# Lists with the variables that have buckets available.
#
rad_buckets = ["ACSWUPT", "ACSWUPTC", "ACSWDNT", "ACSWDNTC", "ACSWUPB",
               "ACSWUPBC", "ACSWDNB", "ACSWDNBC", "ACLWUPT", "ACLWUPTC",
               "ACLWDNT", "ACLWDNTC", "ACLWUPB" ,"ACLWUPBC", "ACLWDNB",
               "ACLWDNBC"]
pr_buckets = ["RAINC", "RAINNC"]


def deaccumulate_var(vardata, varname, wnfiles, wntimes, options):
    #
    # De-accumulates any variable and returns the array of data.
    # By default it deaccumulates it backwards. Uses buckets if
    # it finds the variable in the lists.
    #
    # Check for buckets
    #
    has_buckets = False
    if varname in pr_buckets:
        if wnfiles.current.BUCKET_MM != -1.:
            has_buckets = True
            bucket_size = wnfiles.current.BUCKET_MM
            nbuckets = wnfiles.current.variables["I_" + varname][:]
    elif varname in rad_buckets:
        if wnfiles.current.BUCKET_J != -1.:
            has_buckets = True
            bucket_size = wnfiles.current.BUCKET_J
            nbuckets = wnfiles.current.variables["I_" + varname][:]
    #
    # Read previous files
    #
    if not wnfiles.prv:
        lastval = vardata[0]
        if has_buckets:
            last_nbuckets = nbuckets[0]
            last_nbuckets.shape = (1,) + last_nbuckets.shape
    else:
        lastval = wnfiles.prv.variables[varname][:][-1]
        if has_buckets:
            last_nbuckets = wnfiles.prv.variables["I_" + varname][:][-1]
            last_nbuckets.shape = (1,) + last_nbuckets.shape
    lastval.shape = (1,) + lastval.shape
    shifted_var = np.concatenate([lastval, vardata[:-1]])
    #
    # Deaccumulate
    #
    if has_buckets:
        shifted_nbuckets = np.concatenate([last_nbuckets, nbuckets[:-1]])
        odata = vardata[:] - shifted_var + bucket_size*(nbuckets - shifted_nbuckets)
    else:
        odata = vardata[:] - shifted_var
    return odata


def screenvar_at_2m(varobj, onc, wnfiles, wntimes, options):
    #
    # Works for any variable defined at 2m with no other transformation.
    #
    incvar = wnfiles.current.variables[varobj.varname]
    copyval = np.reshape(incvar[:], incvar.shape[:1]+(1,)+incvar.shape[1:])
    oncvar = get_oncvar(varobj, incvar, onc, options, screenvar_at_2m=True)
    return oncvar, copyval


def screenvar_at_10m(varobj, onc, wnfiles, wntimes, options):
    #
    # Works for any variable defined at 10m with no other transformation.
    #
    incvar = wnfiles.current.variables[varobj.varname]
    copyval = np.reshape(incvar[:],incvar.shape[:1]+(1,)+incvar.shape[1:])
    oncvar = get_oncvar(varobj, incvar, onc, options,
                        screenvar_at_10m=True)
    return oncvar, copyval


def screenvar_at_100m(varobj, onc, wnfiles, wntimes, options):
        #
        # Works for any variable defined at 100m with no other transformation.
        #
        incvar = wnfiles.current.variables[varobj.varname]
        copyval = np.reshape(incvar[:],incvar.shape[:1]+(1,)+incvar.shape[1:])
        oncvar = get_oncvar(varobj, incvar, onc, options,
                            screenvar_at_100m=True)
        return oncvar, copyval


def mask_sea(varobj, onc, wnfiles, wntimes, options):
    #
    # Sets land points to missing.
    #
    incvar = wnfiles.current.variables[varobj.varname]
    if wnfiles.current.variables.has_key("LANDMASK"):
        landmask = wnfiles.current.variables["LANDMASK"][:]
    else:
        if not wnfiles.geo:
            raise IOError(
                "Error: The geo_em file is needed to read the landmask.")
        else:
            landmask = wnfiles.geo.variables["LANDMASK"][:]
    landmask = np.resize(landmask, incvar.shape)
    copyval = np.where(landmask == 1, -9.e+33, incvar[:])
    oncvar = get_oncvar(varobj, incvar, onc, options)
    oncvar.missing_value = np.array(-9.e+33).astype(oncvar.dtype)
    return oncvar, copyval


def mask_land(varobj, onc, wnfiles, wntimes, options):
    #
    # Sets sea points to missing.
    #
    incvar = wnfiles.current.variables[varobj.varname]
    if wnfiles.current.variables.has_key("LANDMASK"):
        landmask = wnfiles.current.variables["LANDMASK"][:]
    else:
        if not wnfiles.geo:
            raise IOError(
                "Error: The geo_em file is needed to read the landmask.")
        else:
            landmask = wnfiles.geo.variables["LANDMASK"][:]
    landmask = np.resize(landmask, incvar.shape)
    copyval = np.where(landmask == 0, -9.e+33, incvar[:])
    oncvar = get_oncvar(varobj, incvar, onc, options)
    oncvar.missing_value = np.array(-9.e+33).astype(oncvar.dtype)
    return oncvar, copyval


def deaccumulate_flux(varobj, onc, wnfiles, wntimes, options):
    #
    # De-accumulates any variable and divides it by the time interval used to
    # deaccumulate in seconds.
    #
    incvar = wnfiles.current.variables[varobj.varname]
    copyval = deaccumulate_var(incvar, varobj.varname, wnfiles, wntimes,
                               options)
    if float(wntimes.outstep_s) == 0:
        raise RuntimeError("Unable to compute a flux with only 1 timestep")

    copyval = copyval/float(wntimes.outstep_s)
    oncvar = get_oncvar(varobj, incvar, onc, options)
    return oncvar, copyval


def deaccumulate(varobj, onc, wnfiles, wntimes, options):
    #
    # De-accumulates any variable if no other transformation is required
    #
    incvar = wnfiles.current.variables[varobj.varname]
    copyval = deaccumulate_var(incvar, varobj.varname, wnfiles, wntimes,
                               options)
    oncvar = get_oncvar(varobj, incvar, onc, options)
    return oncvar, copyval


def fake_extreme(varobj, onc, wnfiles, wntimes, options):
    #
    # Extracts a variable, changes it's name and adds an attribute so extremes
    # must be computed LATER the CDO. It's called "fake extreme" because it is
    # made to replace a extreme that's not computed by CLWRF.
    #
    incvar = wnfiles.current.variables[varobj.varname[:-4]] #[:-4]removes "FMAX" or "FMIN"
    copyval = incvar[:]
    oncvar = get_oncvar(varobj, incvar, onc, options)
    oncvar.warning = ("This is not a real extreme, extremes still need to be "
                      "computed.")
    return oncvar, copyval


def fake_extreme_screenvar_at_2m(varobj, onc, wnfiles, wntimes, options):
    #
    # Extracts a variable, changes it's name and adds an attribute so extremes
    # must be computed LATER the CDO. It's called "fake extreme" because it is
    # made to replace a extreme that's not computed by CLWRF.
    #
    incvar = wnfiles.current.variables[varobj.varname[:-4]] #[:-4]removes "MAX" or "MIN"
    copyval = np.reshape(incvar[:],incvar.shape[:1]+(1,)+incvar.shape[1:])
    oncvar = get_oncvar(varobj, incvar, onc, options, screenvar_at_2m=True)
    oncvar.warning = ("This is not a real extreme, extremes still need to be "
                      "computed.")
    return oncvar, copyval


def fake_extreme_screenvar_at_10m(varobj, onc, wnfiles, wntimes, options):
    #
    # Extracts a variable, changes it's name and adds an attribute so extremes
    # must be computed LATER the CDO. It's called "fake extreme" because it is
    # made to replace a extreme that's not computed by CLWRF.
    #
    incvar = wnfiles.current.variables[varobj.varname[:-4]] #[:-4]removes "MAX" or "MIN"
    copyval = np.reshape(incvar[:],incvar.shape[:1]+(1,)+incvar.shape[1:])
    oncvar = get_oncvar(varobj, incvar, onc, options, screenvar_at_10m=True)
    oncvar.warning = ("This is not a real extreme, extremes still need to be "
                      "computed.")
    return oncvar, copyval


def fake_extreme_screenvar_at_100m(varobj, onc, wnfiles, wntimes, options):
        #
        # Extracts a variable, changes it's name and adds an attribute so
        # extremes must be computed LATER the CDO. It's called "fake extreme"
        # because it is made to replace a extreme that's not computed by CLWRF.
        #
        incvar = wnfiles.current.variables[varobj.varname[:-4]] #[:-4]removes "MAX" or "MIN"
        copyval = np.reshape(incvar[:],incvar.shape[:1]+(1,)+incvar.shape[1:])
        oncvar = get_oncvar(varobj, incvar, onc, options,
                            screenvar_at_100m=True)
        oncvar.warning = "This is not a real extreme, extremes still need to be computed."
        return oncvar, copyval


def rotate_uas(varobj, onc, wnfiles, wntimes, options):
    uvarname = varobj.varname[:-2] # remove the "ER"
    vvarname = "V" + varobj.varname[1:-2]  # remove the "U" and the "ER"
    u = wnfiles.current.variables[uvarname]
    v = wnfiles.current.variables[vvarname]
    if not wnfiles.geo:
        raise IOError("Error: The geo_em file is needed to rotate the winds")
    else:
        sina = wnfiles.geo.variables["SINALPHA"][:]
        cosa = wnfiles.geo.variables["COSALPHA"][:]
    copyval = (
        u[:]*cosa[np.newaxis,...] - v[:]*sina[np.newaxis,...]
    )
    copyval.shape = u.shape[:1]+ (1,) + u.shape[1:]
    oncvar = get_oncvar(varobj, u, onc, options, screenvar_at_10m=True)
    return oncvar, copyval


def rotate_vas(varobj, onc, wnfiles, wntimes, options):
    uvarname = "U" + varobj.varname[1:-2] # remove the V and the "ER"
    vvarname = varobj.varname[:-2] # remove the "ER"
    u = wnfiles.current.variables[uvarname]
    v = wnfiles.current.variables[vvarname]
    if not wnfiles.geo:
        raise IOError("Error: The geo_em file is needed to rotate the winds")
    else:
        sina = wnfiles.geo.variables["SINALPHA"][:]
        cosa = wnfiles.geo.variables["COSALPHA"][:]
    copyval = (
        u[:]*sina[np.newaxis,...] + v[:]*cosa[np.newaxis,...]
    )
    copyval.shape = v.shape[:1]+ (1,) + v.shape[1:]
    oncvar = get_oncvar(varobj, v, onc, options, screenvar_at_10m=True)
    return oncvar, copyval


def rotate_uu(varobj, onc, wnfiles, wntimes, options):
        uvarname = "UU"
        vvarname = "VV"
        u = wnfiles.current.variables[uvarname]
        v = wnfiles.current.variables[vvarname]
        if not wnfiles.geo:
                raise IOError(
                    "Error: The geo_em file is needed to rotate the winds")
        else:
                sina = wnfiles.geo.variables["SINALPHA"][:]
                cosa = wnfiles.geo.variables["COSALPHA"][:]
        copyval = (
                u[:]*cosa[np.newaxis,...] - v[:]*sina[np.newaxis,...]
        )
        copyval.shape = u.shape[:1]+ (1,) + u.shape[1:]
        oncvar = get_oncvar(varobj, u, onc, options, screenvar_at_100m=True)
        return oncvar, copyval


def rotate_vv(varobj, onc, wnfiles, wntimes, options):
        uvarname = "UU"
        vvarname = "VV"
        u = wnfiles.current.variables[uvarname]
        v = wnfiles.current.variables[vvarname]
        if not wnfiles.geo:
                raise IOError(
                    "Error: The geo_em file is needed to rotate the winds")
        else:
                sina = wnfiles.geo.variables["SINALPHA"][:]
                cosa = wnfiles.geo.variables["COSALPHA"][:]
        copyval = (
                u[:]*sina[np.newaxis, ...] + v[:]*cosa[np.newaxis, ...]
        )
        copyval.shape = v.shape[:1] + (1,) + v.shape[1:]
        oncvar = get_oncvar(varobj, v, onc, options, screenvar_at_100m=True)
        return oncvar, copyval