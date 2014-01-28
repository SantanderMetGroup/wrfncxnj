from wrfncxnj_base import *
from wrfncxnj_fun import *

if opt.quiet:
	print "Extract and join is running in quiet mode"
wrfncIter = WrfNcFiles()
#
# Apply the options that have been parsed.
#
if opt.geofile:
	wrfncIter.geo = ncdf.Dataset(opt.geofile, "r")
if opt.fullfile:
	wrfncIter.full = ncdf.Dataset(opt.fullfile, "r")

if not opt.OFILE and not opt.OUTPUT_PATTERN:
	sys.stderr.write("Missing output file or pattern!")
	sys.exit(1)
if opt.OFILE and opt.OUTPUT_PATTERN:
	sys.stderr.write("-o and --output-pattern cannot be used at the same time!")
	sys.exit(1)
if opt.splitvars and not opt.OUTPUT_PATTERN:
	sys.stderr.write("Split variables has been activated, but an output pattern has not been specified. Use --output-pattern= to specify one")
	sys.exit(1)
if not opt.splitvars and opt.OUTPUT_PATTERN:
	sys.stderr.write("An output pattern is needed to use --split-in-variables. Please set it with the --output-pattern option")
	sys.exit(1)
if not opt.splitvars and opt.selected_plevs:
	sys.stderr.write("--split-in-vars and -sel-plevs must be used at the same time")
	sys.exit(1)
ofh = OfileHandler()
ofh.set_tempdir()

if opt.globfiles:
	files = glob(opt.globfiles)
	files.sort()
elif opt.filelist:
	files = map(string.strip, open(opt.filelist, "r").readlines())
#elif opt.singlerec:
#	files = [opt.geofile]
else:
	files = args
vars = opt.vars.split(',')
if not opt.vtable:
	opt.vtable = sys.argv[0].replace(".py", ".table")
#
# Look for geofile
#
is_geofile = False
if not files and wrfncIter.geo:
	print "No input files provided."
	print "Trying to find the variables in the geo_em file provided."
	files = [opt.geofile,]
	is_geofile = True
#
# Get the projection
#
proj = Projection()

proj.get_projection(files[0])
#
# Define the list of variable objects.
# 
vars = stdvars(vars, opt.vtable, proj, files[0])
#
#	Clone the structure of the netcdf file and get the initial time from the first file.
#
if not opt.refdate:
	opt.refdate = "1950-01-01_00:00:00"
wnt = WrfNcTime(datetime.strptime(opt.refdate, '%Y-%m-%d_%H:%M:%S'))
#
# Create output netCDF files
#
oncs, oncnames = create_oncs(vars, ofh, files, wnt, proj)
ofh.set_oncnames(oncnames)
#
# Get coordinate variable objects from the output netCDF files.
#
onclon, onclat, onctime	= get_onc_lon_lat_time(oncs)
#
# Loop over files extracting variables and times
#
if opt.discard:
	files = discard_suspect_files(files, opt.discard)
wrfncIter.loadFiles(files, opt.prevfile, opt.nextfile)
for wrfnc in wrfncIter:
	if not opt.quiet: print wrfnc
	if not opt.singlerec:
		wnt.checkStep(wrfnc)
	times = wnt.getTimes(wrfnc, opt.singlerec)
	if not opt.quiet: print wnt
	#
	#	Set times and loop variables
	#
	for varname in vars:
		if not opt.splitvars:
			onc = oncs
		else:
			onc = oncs[varname]
		if not opt.quiet: print "Processing var %s" % varname 

		if "compute_%s" % varname in locals():
			if not opt.quiet: print "	Using the compute_%s function" % varname
			process_func = locals()["compute_%s" % varname]
			oncvar, copyval = process_func(vars[varname], onc, wrfnc, wnt)
			oncvar.set_ncvars(vars[varname])
			oncvar.assign_values(copyval[:wnt.nrec].astype('f'), wnt.iini, wnt.iend, times)
			oncvar.sync()		
		elif vars[varname].transform:
			pt = ParseTransform(vars[varname].transform)
			if not opt.quiet: print pt
			vardic = {}
			for var in pt.variables:
				vardic[var] = wrfnc.current.variables[var][:]
			oncvar, copyval = pt.execute(vars[varname], onc, wrfnc, wnt, vardic)
			oncvar.set_ncvars(vars[varname])
			oncvar.assign_values(copyval[:wnt.nrec].astype('f'), wnt.iini, wnt.iend, times)
			oncvar.sync()
		else:
			incvar = wrfnc.current.variables[varname]
			oncvar = get_oncvar(vars[varname], incvar, onc)
			oncvar.set_ncvars(vars[varname])
			oncvar.assign_values(incvar[:wnt.nrec].astype('f'), wnt.iini, wnt.iend, times)
			oncvar.sync()
	wnt.cycle()

oncvar.sync()
oncvar.close()
if opt.ftimes:
  time1, time2 = opt.ftimes.split(",")
  ofh.filter_times(time1, time2)
if opt.TEMPDIR:
		ofh.move_output()
wrfncIter.closeCommon()
print "Extract and join finished successfully"
