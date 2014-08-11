# -*- coding: utf-8 -*-
import numpy as np
from wrfncxnj_base import get_oncvar, Constants, opt, compute_mslp, deaccumulate_var
import sys

def compute_temperature(p, t):
	# Some required physical constants:
	return (t+300.)*(p/Constants.p1000mb)**Constants.rcp

def compute_CLTFR(varobj, onc, wnfiles, wntimes):
	"""Computation of total cloud fraction in case CLT	(from p_interp) is not available
	cldfra: partial cloud fraction [1]

	Following Sundqvist, 1989, Mont. Weather Rev.
	NOTE: This computation is not really comparable to CLT, since is
	computed from pressure levels instead of the native sigma levels
	"""
	if wnfiles.current.variables.has_key("CLT"): # The file was correctly processed by p_interp
		var = wnfiles.current.variables["CLT"]
		copyval = var[:]
		oncvar = get_oncvar(varobj, var, onc)
	else:
		# We need to compute it using Sundqvist (1989) MWR, eq 5.1
		if not opt.quiet: print "Computing total cloud cover in Python!"
		cldfra =	wnfiles.current.variables["CLDFRA"]
		bj = cldfra[:]
		ceros = np.zeros(cldfra.shape[:1]+(1,)+cldfra.shape[2:],cldfra.dtype)
		bjm1 = np.concatenate([ceros, cldfra[:,:-1,:,:]], axis=1)
		bjmax = np.maximum(bjm1,bj)
		denom = 1.-bjm1
		denom = np.where(denom == 0., 1., denom) # Avoid singularities
		copyval = 1. - np.multiply.reduce((1.-bjmax)/denom, axis=1)
		oncvar = get_oncvar(varobj, cldfra, onc, out_is_2D_but_in_3D=True)
	return oncvar, copyval

def compute_RH2(varobj, onc, wnfiles, wntimes):
	"""Computation of relative humidity at 2 m 
	t2: 2m temperature [K]
	q2: 2m mixing ratio [kg kg-1]
	psfc: surface prssure (assumed the same at 2m) [Pa]
	e: vapor pressure in air [Pa]
	es: saturation vapor pressure [Pa]

	Some required physical constants: epsilon_gamma, tkelvin, es_[A/B]tetens_[ice/vapor]
	"""
	t2 =	wnfiles.current.variables["T2"][:]
	q2 =	wnfiles.current.variables["Q2"][:]
	psfc =	wnfiles.current.variables["PSFC"]
	e = q2*psfc[:]/(100.*(Constants.epsilon_gamma+q2))
	es = np.where(
		t2-Constants.tkelvin <= 0., 
		Constants.es_base_tetens*10.**(((t2-Constants.tkelvin)*Constants.es_Atetens_ice)/
			((t2-Constants.tkelvin)+Constants.es_Btetens_ice)), 
		Constants.es_base_tetens*10.**(((t2-Constants.tkelvin)*Constants.es_Atetens_vapor)/
			((t2-Constants.tkelvin)+Constants.es_Btetens_vapor)))
	copyval = e/es
	copyval.shape=psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, psfc, onc, screenvar_at_2m=True)
#	rh2 = e/es
#	copyval = max(min(rh2, 1.), 0.)
	return oncvar, copyval

def compute_RH2MEAN(varobj, onc, wnfiles, wntimes):
	"""Computation of relative humidity at 2 m 
	t2: 2m temperature [K]
	q2: 2m mixing ratio [kg kg-1]
	psfc: surface prssure (assumed the same at 2m) [Pa]
	e: vapor pressure in air [Pa]
	es: saturation vapor pressure [Pa]

	Some required physical constants: epsilon_gamma, tkelvin, es_[A/B]tetens_[ice/vapor]
	"""
	t2 =	wnfiles.current.variables["T2MEAN"][:]
	q2 =	wnfiles.current.variables["Q2MEAN"][:]
	psfc =	wnfiles.current.variables["PSFCMEAN"]
	e = q2*psfc[:]/(100.*(Constants.epsilon_gamma+q2))
	es = np.where(
		t2-Constants.tkelvin <= 0.,
		Constants.es_base_tetens*10.**(((t2-Constants.tkelvin)*Constants.es_Atetens_ice)/
			((t2-Constants.tkelvin)+Constants.es_Btetens_ice)),
		Constants.es_base_tetens*10.**(((t2-Constants.tkelvin)*Constants.es_Atetens_vapor)/
			((t2-Constants.tkelvin)+Constants.es_Btetens_vapor)))
	copyval = e/es
	copyval.shape=psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, psfc, onc, screenvar_at_2m=True)
#	rh2 = e/es
#	copyval = max(min(rh2, 1.), 0.)
	return oncvar, copyval

def compute_SPDUV10(varobj, onc, wnfiles, wntimes):
	"""Computation of wind speed at 10 m
	u10: 10 m horizontal WE-wind speed
	v10: 10 m horizontal NS-wind speed
	"""
	u10 =	wnfiles.current.variables["U10"]
	v10 =	wnfiles.current.variables["V10"][:]
	copyval = np.sqrt(u10[:]*u10[:] + v10*v10)
	copyval.shape=u10.shape[:1]+(1,)+u10.shape[1:]
	oncvar = get_oncvar(varobj, u10, onc, screenvar_at_2m=True)
	return oncvar, copyval

def compute_TDPS(varobj, onc, wnfiles, wntimes):
	"""Computation of relative humidity at 2 m 
	t2: 2m temperature [K]
	q2: 2m mixing ratio [kg kg-1]
	psfc: surface prssure (assumed the same at 2m) [Pa]
	e: vapor pressure in air [hPa]
	es: saturation vapor pressure [hPa]

	Some required physical constants: tkelvin, es_base_bolton, es_[A/B]bolton, RdRv
	"""
	t2 =	wnfiles.current.variables["T2"][:]
	q2 =	wnfiles.current.variables["Q2"][:]
	psfc =	wnfiles.current.variables["PSFC"]
	es = 10.* Constants.es_base_bolton * np.exp(Constants.es_Abolton*(t2-Constants.tkelvin)
		/(t2-Constants.tkelvin+Constants.es_Bbolton))
	rh = 0.01*psfc[:]/es*(q2/(Constants.RdRv+q2))
	rh = np.where(rh > 1., 1., rh)
	rh = np.where(rh < 0., 0., rh)

	e = rh*Constants.es_base_bolton*np.exp(Constants.es_Abolton*(t2-Constants.tkelvin)
		/(t2-Constants.tkelvin+Constants.es_Bbolton))
	copyval = (116.9+237.3*np.log(e))/(16.78-np.log(e))+Constants.tkelvin
	copyval.shape=psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, psfc, onc, screenvar_at_2m=True)
	return oncvar, copyval

def compute_TDPS_its90(varobj, onc, wnfiles, wntimes):
	"""Computation of dew point temperature at 2 m
	q2: 2m mixing ratio [kg kg-1]
	psfc: surface prssure (assumed the same at 2m) [Pa]
	lnes: pure saturated vapor pressure [Pa]
	Following ITS-90 formulation

	Some required physical constants: epsilon_gamma, es_[A/B]tetens_vapor
	"""
	t2 =	wnfiles.current.variables["T2"][:]
	psfc =	wnfiles.current.variables["PSFC"]

	g_es0 = -2.8365744e3
	g_es1 = -6.028076559e3
	g_es2 = 1.954263612e1
	g_es3 = -2.737830188e-2
	g_es4 = 1.6261698e-5
	g_es5 = 7.0229056e-10
	g_es6 = -1.8680009e-13
	g_es7 = 2.7150305

	k_es0 = -5.8666426e3
	k_es1 = 2.232870244e1
	k_es2 = 1.39387003e-2
	k_es3 = -3.4262402e-5
	k_es4 = 2.7040955e-8
	k_es5 = 6.7063522e-1

	c_tdv0 = 2.0798233e2
	c_tdv1 = -2.0156028e1
	c_tdv2 = 4.6778925e-1
	c_tdv3 = -9.2288067e-6
	d_tdv0 = 1.
	d_tdv1 = -1.3319669e-1
	d_tdv2 = 5.6577518e-3
	d_tdv3 = -7.5172865e-5

	c_tdi0 = 2.1257969e2
	c_tdi1 = -1.0264612e1
	c_tdi2 = 1.4354796e-1
	d_tdi0 = 1.
	d_tdi1 = -8.2871619e-2
	d_tdi2 = 2.3540411e-3
	d_tdi3 = -2.4363951e-5

	lnes = np.where(
		t2-Constants.tkelvin <= 0.,
		k_es0*t2**(-1)+ k_es1*t2**(0)+ k_es2*t2**(1)+ k_es3*t2**(2)+ k_es4*t2**(3)+ 
		k_es5*np.log(t2),
		g_es0*t2**(-2)+ g_es1*t2**(-1)+ g_es2*t2**(0)+ g_es3*t2**(1)+ g_es4*t2**(2)+		
		g_es5*t2**(3)+ g_es6*t2**(4)+ g_es7*np.log(t2))
	td2A = np.where(
		t2-Constants.tkelvin <= 0.,
		c_tdi0*lnes**0+ c_tdi1*lnes**1+ c_tdi2*lnes**2,
		c_tdv0*lnes**0+ c_tdv1*lnes**1+ c_tdv2*lnes**2+ c_tdv3*lnes**3)
	td2B = np.where(
		t2-Constants.tkelvin <= 0.,
		d_tdi0*lnes**0+ d_tdi1*lnes**1+ d_tdi2*lnes**2+ d_tdi3*lnes**3,
		d_tdv0*lnes**0+ d_tdv1*lnes**1+ d_tdv2*lnes**2+ d_tdv3*lnes**3)
 
	copyval = td2A/td2B
	copyval.shape=psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, psfc, onc, screenvar_at_2m=True)
	return oncvar, copyval

def compute_UER(varobj, onc, wnfiles, wntimes):
	if wnfiles.current.variables.has_key("UU"):	# wind on p-levels from p_interp
		u = wnfiles.current.variables["UU"]
		v = wnfiles.current.variables["VV"]
	else:
		u = wnfiles.current.variables["U"]
		v = wnfiles.current.variables["V"]
	if not wnfiles.geo:
		print >> sys.stderr, "Error: The geo_em file is needed to rotate the winds"
		sys.exit(1)
	else:
		sina = wnfiles.geo.variables["SINALPHA"][:]
		cosa = wnfiles.geo.variables["COSALPHA"][:]
	copyval = (
		u[:]*cosa[np.newaxis,...] - v[:]*sina[np.newaxis,...]
	)
	oncvar = get_oncvar(varobj, u, onc)
	return oncvar, copyval

def compute_VER(varobj, onc, wnfiles, wntimes):
	if wnfiles.current.variables.has_key("UU"):	# wind on p-levels from p_interp
		u = wnfiles.current.variables["UU"]
		v = wnfiles.current.variables["VV"]
	else:
		u = wnfiles.current.variables["U"]
		v = wnfiles.current.variables["V"]
	if not wnfiles.geo:
		print >> sys.stderr, "Error: The geo_em file is needed to rotate the winds"
		sys.exit(1)
	else:
		sina = wnfiles.geo.variables["SINALPHA"][:]
		cosa = wnfiles.geo.variables["COSALPHA"][:]
	copyval = (
		u[:]*sina[np.newaxis,...] + v[:]*cosa[np.newaxis,...]
	)
	oncvar = get_oncvar(varobj, v, onc)
	return oncvar, copyval



def compute_SMOIS1(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables['SMOIS']
	if wnfiles.current.variables.has_key("DZS"):
		layer_width = wnfiles.current.variables['DZS'][:]
	else:
		if wnfiles.full.variables.has_key("DZS"):
			layer_width = wnfiles.full.variables['DZS'][:]
		else:	
			print >> sys.stderr, "Error: 'DZS' variable not found!"
			sys.exit(1)
	smois1 = incvar[:,0,:,:]*layer_width[0,0]*1000. # m3/m3 -> Kg/m2
	# Sets sea points to missing values.
	if wnfiles.current.variables.has_key("LANDMASK"):
		landmask = wnfiles.current.variables["LANDMASK"][:]
	else:
		if not wnfiles.geo:
			print >> sys.stderr, "Error: The geo_em file is needed to read the landmask."
			sys.exit(1)
		else:
			landmask = wnfiles.geo.variables["LANDMASK"][:]
	landmask = np.resize(landmask, smois1.shape)
	copyval = np.where( landmask == 0, -9.e+33, smois1[:]) 
	oncvar = get_oncvar(varobj, incvar, onc, out_is_2D_but_in_3D=True)
	oncvar.missing_value = np.array(-9.e+33).astype(oncvar.dtype)
	return oncvar, copyval

def compute_MRSO(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables['SMOIS']
	if wnfiles.current.variables.has_key("DZS"):
		layer_width = wnfiles.current.variables['DZS'][0]
	else:
		if wnfiles.full.variables.has_key("DZS"):
			layer_width = wnfiles.full.variables['DZS'][0]
		else:
			print >> sys.stderr, "Error: 'DZS' variable not found!"
			sys.exit(1)
	smois = incvar[:]
	smois_byarea = smois*layer_width[np.newaxis,:,np.newaxis,np.newaxis]		 
	smois_total = np.sum(smois_byarea, axis=1)*1000.
	if wnfiles.current.variables.has_key("LANDMASK"):
		landmask = wnfiles.current.variables["LANDMASK"][:]
	else:
		if not wnfiles.geo:
			print  >> sys.stderr, "Error: I need the geo_em file to read the landmask!"
			sys.exit(1)
		else:
			landmask = wnfiles.geo.variables["LANDMASK"][:]
	landmask = np.resize(landmask, smois_total.shape)
	copyval = np.where( landmask == 0, -9.e+33, smois_total[:]) 
	oncvar = get_oncvar(varobj, incvar, onc, out_is_2D_but_in_3D=True)
	oncvar.missing_value = np.array(-9.e+33).astype(oncvar.dtype)
	return oncvar, copyval

def compute_RAINF(varobj, onc, wnfiles, wntimes):
	"""Deaccumulates the precipitation field

	This function adds RAINNC and RAINC, support for RAINTOT has been deprecated since 
	it does not support buckets. It deaccumulates from the value on the previous output
	time step. The flux is computed dividing by the time interval in seconds.
	"""
	incvar = wnfiles.current.variables["RAINNC"]
	rainnc = incvar[:]
	rainc = wnfiles.current.variables["RAINC"][:]
	deac_rainnc = deaccumulate_var(rainnc, "RAINNC", wnfiles, wntimes)
	deac_rainc  = deaccumulate_var(rainc, "RAINC", wnfiles, wntimes)

	copyval = deac_rainnc + deac_rainc
	#copyval = np.where(copyval<0., 0, copyval)/float(wntimes.outstep_s)
	copyval = copyval/float(wntimes.outstep_s)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_MRRO(varobj, onc, wnfiles, wntimes):
        """Deaccumulates the total runoff field

        This function adds SFROFF and UDROFF. It deaccumulates from the value on the previous 
        output time step. The flux is computed dividing by the time interval in seconds.
        """
        incvar = wnfiles.current.variables["SFROFF"]
        sfroff = incvar[:]
        udroff = wnfiles.current.variables["UDROFF"][:]
        deac_sfroff = deaccumulate_var(sfroff, "SFROFF", wnfiles, wntimes)
        deac_udroff  = deaccumulate_var(udroff, "UDROFF", wnfiles, wntimes)

        copyval = deac_sfroff + deac_udroff
        #copyval = np.where(copyval<0., 0, copyval)/float(wntimes.outstep_s)
        copyval = copyval/float(wntimes.outstep_s)
        oncvar = get_oncvar(varobj, incvar, onc)
        return oncvar, copyval

def compute_RAIN(varobj, onc, wnfiles, wntimes):
	"""Deaccumulates the precipitation field
	This function adds RAINNC and RAINC, RAINTOT has been deprecated since it does not support buckets.
 It deaccumulates from the value on the previous output time step.
	"""	
	incvar = wnfiles.current.variables["RAINNC"]
	rainnc = incvar[:]
	rainc = wnfiles.current.variables["RAINC"][:]
	deac_rainnc = deaccumulate_var(rainnc, "RAINNC", wnfiles, wntimes)
	deac_rainc  = deaccumulate_var(rainc, "RAINC", wnfiles, wntimes)

	copyval = deac_rainnc + deac_rainc
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_RAINFORWARD(varobj, onc, wnfiles, wntimes):
	if wnfiles.current.variables.has_key("RAINTOT"): # The file was processed by p_interp
		incvar = wnfiles.current.variables["RAINTOT"]
		pr = incvar[:]
	else:	# We should add convective and large-scale rainfall
		incvar = wnfiles.current.variables["RAINNC"]
		pr = incvar[:] + wnfiles.current.variables["RAINC"][:]
	if not wnfiles.nxt:
		nextpr = pr[-1]
	elif wnfiles.current.variables.has_key("RAINTOT"):
		nextpr = wnfiles.nxt.variables["RAINTOT"][0]
	else:
		nextpr = wnfiles.nxt.variables["RAINNC"][0] + wnfiles.nxt.variables["RAINC"][0]
	nextpr.shape = (1,) + nextpr.shape
	copyval = np.concatenate([pr[1:], nextpr]) - pr
	copyval = np.where(copyval<0., 0, copyval)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_ACRLS(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["ACLWUPB"]
	aclwupb = incvar[:]
	deac_aclwupb = deaccumulate_var(aclwupb, "ACLWUPB", wnfiles, wntimes)
	#
	# If there is not ACLWDNB in the file it deaccumulates just ACLWUPB and the uses GLW to compute ACRLS.
	#
	if wnfiles.current.variables.has_key("ACLWDNB"):
		aclwdnb = wnfiles.current.variables["ACLWDNB"][:]
		deac_aclwdnb = deaccumulate_var(aclwdnb, "ACLWDNB", wnfiles, wntimes)
		deac_rlt =	deac_aclwdnb - deac_aclwupb
		copyval = deac_rlt/float(wntimes.outstep_s)
	else:
		glw = wnfiles.current.variables["GLW"][:]
		if not wnfiles.prv:
			glw[0] = 0 #If there is not previous file makes 0 the first glw, so the substraction below is also 0.
		copyval = (glw - deac_aclwupb)/float(wntimes.outstep_s)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_RST(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["ACSWUPT"]
	acswupt = incvar[:]
	deac_acswupt = deaccumulate_var(acswupt, "ACSWUPT", wnfiles, wntimes)
	acswdnt = wnfiles.current.variables["ACSWDNT"][:]
	deac_acswdnt = deaccumulate_var(acswdnt, "ACSWDNT", wnfiles, wntimes)
	rst = deac_acswdnt - deac_acswupt
	copyval = rst/float(wntimes.outstep_s)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_RLB(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["ACLWUPB"]
	aclwupb = incvar
	deac_aclwupb = deaccumulate_var(aclwupb, "ACLWUPB", wnfiles, wntimes)
	aclwdnb = wnfiles.current.variables["ACLWDNB"]
	deac_aclwdnb = deaccumulate_var(aclwdnb, "ACLWDNB", wnfiles, wntimes)

	copyval = (deac_aclwupb - deac_aclwdnb)/float(wntimes.outstep_s)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_RLT(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["OLR"]
	olr = incvar
	aclwdnt = wnfiles.current.variables["ACLWDNT"]
	deac_aclwdnt = deaccumulate_var(aclwdnt, "ACLWDNT", wnfiles, wntimes)
	copyval = (olr - deac_aclwdnt)/float(wntimes.outstep_s)
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_WINDFMAX(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["U10"]
	u = incvar[:]
	v = wnfiles.current.variables["V10"][:]
	copyval = np.sqrt(u*u+v*v)
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_10m=True)
	oncvar.warning = "This is not a real extreme, extremes still need to be computed." 
	return oncvar, copyval

def compute_RAINFMAX1H(varobj, onc, wnfiles, wntimes):
	"""Deaccumulates the precipitation field and changes it's name 

	This function looks for RAINTOT if available, otherwise adds up RAINNC and
	RAINCV. It deaccumulates from the value on the previous output time step.
	A flux is computed dividing by the timestep in seconds.
	"""
	if wnfiles.current.variables.has_key("RAINTOT"): # The file was processed by p_interp
		incvar = wnfiles.current.variables["RAINTOT"]
		pr = incvar[:]
	else:	# We should add convective and large-scale rainfall
		incvar = wnfiles.current.variables["RAINNC"]
		pr = incvar[:] + wnfiles.current.variables["RAINC"][:]
	if not wnfiles.prv:
		lastpr = pr[0]
	elif wnfiles.current.variables.has_key("RAINTOT"):
		lastpr = wnfiles.prv.variables["RAINTOT"][:][-1]
	else:
		lastpr = wnfiles.prv.variables["RAINNC"][:][-1] + wnfiles.prv.variables["RAINC"][:][-1]
	lastpr.shape = (1,) + lastpr.shape
	copyval = pr - np.concatenate([lastpr,pr[:-1]])
	copyval = np.where(copyval<0., 0, copyval)
	oncvar = get_oncvar(varobj, incvar, onc)
	oncvar.warning = "This is not a real extreme, extremes still need to be computed."
	return oncvar, copyval
def compute_RHODRY2(varobj, onc, wnfiles, wntimes):
	"""Computes the 2 meter dry air density from 2m temperature, and surface pressure.
	"""
	incvar = wnfiles.current.variables["T2"] # K
	tas = incvar[:]
	psfc = wnfiles.current.variables["PSFC"][:] # Pa
	Rd = 286.9 # J kg -1 K -1
	rho = (psfc/(Rd*tas))[:]
	copyval = rho[:]
	copyval.shape = psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_2m=True)
	return oncvar, copyval
def compute_RHO2(varobj, onc, wnfiles, wntimes):
	"""Computes the 2 meter dry air density from 2m temperature, and surface pressure.
	"""
	incvar = wnfiles.current.variables["T2"] # K
	tas = incvar[:]
	if wnfiles.current.variables.has_key("PSFC"):
		psfc = wnfiles.current.variables["PSFC"][:] # Pa
	else:
		mslp = wnfiles.current.variables["MSLP"][:] # Pa
		height = wnfiles.geo.variables["HGT_M"][:] # m
		psfc = mslp - height*10 #10 Pa/m gradient

	q = wnfiles.current.variables["Q2"][:] # kg kg -1
	Rd = 286.9 # J kg -1 K -1
	tasv = (tas*(1 + 0.61*q))[:]
	rho = (psfc/(Rd*tasv))[:]
	copyval = rho[:]
	copyval.shape = psfc.shape[:1]+(1,)+psfc.shape[1:]
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_2m=True)
	return oncvar, copyval
# old wrfncxnj.py elif blocks variables
def compute_WIND10(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["U10"]
	u = incvar[:]
	v = wnfiles.current.variables["V10"][:]
	copyval = np.sqrt(u*u + v*v)
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_10m=True)
	return oncvar, copyval

def compute_U10(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["U10"]
	copyval = np.reshape(incvar[:wntimes.nrec], incvar.shape[:1]+(1,)+incvar.shape[1:])
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_10m=True)
	return oncvar, copyval

def compute_V10(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables["V10"]
	copyval = np.reshape(incvar[:wntimes.nrec], incvar.shape[:1]+(1,)+incvar.shape[1:])
	oncvar = get_oncvar(varobj, incvar, onc, screenvar_at_10m=True)
	return oncvar, copyval

def compute_GEOP(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables['PH']
	copyval = incvar[:wntimes.nrec] + wnfiles.current.variables["PHB"][:wntimes.nrec]
	# De-stagger the geopotential
	copyval = (copyval[:,:-1]+copyval[:,1:])/2.
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_TEMP(varobj, onc, wnfiles, wntimes):
	if wnfiles.current.variables.has_key("TT"):
		incvar = wnfiles.current.variables['TT']
		copyval = incvar
	else:
		incvar = wnfiles.current.variables['T']
		pres = wnfiles.current.variables['P'][:wntimes.nrec] + wnfiles.current.variables["PB"][:wnt.nrec]
		copyval = compute_temperature(pres, incvar[:wntimes.nrec])
	oncvar = get_oncvar(varobj, incvar, onc)
	return oncvar, copyval

def compute_MSLP(varobj, onc, wnfiles, wntimes):
	if wnfiles.current.variables.has_key("MSLP"):
		incvar = wnfiles.current.variables['MSLP']
		copyval = incvar[:wntimes.nrec]
		oncvar = get_oncvar(varobj, incvar, onc)
	else:
		incvar = wnfiles.current.variables['P']
		p = incvar[:]
		pb = wnfiles.current.variables['PB'][:]
		ph = wnfiles.current.variables['PH'][:]
		phb = wnfiles.current.variables['PHB'][:]
		t = wnfiles.current.variables['T'][:]
		qvapor = wnfiles.current.variables['QVAPOR'][:]
		copyval = compute_mslp(p, pb, ph, phb, t , qvapor)
		oncvar = get_oncvar(varobj, incvar, onc, out_is_2D_but_in_3D=True)
        return oncvar, copyval

def compute_VIS(varobj, onc, wnfiles, wntimes):
	incvar = wnfiles.current.variables['QCLOUD']
	qcloud = incvar[:, 0, :, :] # kg kg-1
	density = 1/wnfiles.current.variables['ALT'][:, 0, :, :] # kg m-3
        qcdensity = 1000*qcloud*density # g m-3
        qcdensity = np.where(qcdensity < 0, 0, qcdensity)
	#
	# Visibility following Stoelinga and Warner (1998)
	#
	beta = 144.7*(qcdensity**0.88) # km-1
        beta = np.where(beta <= 1e-20, 1e-20, beta)
        vis = np.empty(beta.shape)
        vis = np.where(beta == 0, 11000, (-np.log(0.02)/beta)*1000) # m 
	#vis = np.where(beta == 0, 11000, vis)
        vis = np.where(vis > 11000, 11000, vis)
	copyval = vis
	oncvar = get_oncvar(varobj, incvar, onc, out_is_2D_but_in_3D=True)
	return oncvar, copyval
