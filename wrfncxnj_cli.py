from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(quiet=False,singlerec=False)
parser.add_option(
	"-f", "--files", dest="globfiles",
	help="Regular expression to be parsed by python to get the input files to process", metavar="REGEXP"
)
parser.add_option(
	"--from-file", dest="filelist", default="",
	help="Text file containing the input files. One per row", metavar="FILELIST.txt"
)
parser.add_option(
	"--previous-file", dest="prevfile", default="",
	help="Extra input file to be prepended to the input files ONLY for BACKWARD deaccumulations", metavar="WRFNCFILE.nc"
)
parser.add_option(
	"--next-file", dest="nextfile", default="",
	help="Extra input file to be prepended to the input files ONLY for FORWARD deaccumulations", metavar="WRFNCFILE.nc"
)
parser.add_option(
	"-v", "--variables", dest="vars",
	help="Variables to extract. Apart from those defined in the file, you can ask for any of the following derived variables: MSLP, U10ER, V10ER, WIND", metavar="VAR1[,VAR2,...]"
)
parser.add_option(
	"-d", "--discard-criteria", dest="discard",
	help="Enable discarding files. Currently only the uncommon_size criteria is implemented", metavar="uncommon_size"
)
parser.add_option(
	"-t", "--variable-table", dest="vtable",
	help="Table for translating WRF names into IPCC standard names", metavar="variable.table"
)
parser.add_option(
	"-a", "--attributes", dest="attributes",
	help="Table for setting the global attributes of the file", metavar="atributes.file"
)
parser.add_option(
	"-q", "--quiet", action="store_true", default = False,
	help="Run quietly"
)
parser.add_option(
	"--single-record", action="store_true", dest="singlerec",
	help="Save only one record. Useful to extract fixed fields (LANDMASK, LANDUSE, ...)"
)
parser.add_option(
	"-z", action="store_true", default=False, dest="zaxis",
	help="Create Z axis information"
)
parser.add_option(
	"-s", action="store_true", default=False, dest="saxis",
	help="Create soil layer axis information"
)
parser.add_option(
	"-p", action="store_true", default=False, dest="paxis",
	help="Create pressure level axis information"
)
parser.add_option(
	"-m", action="store_true", default=False, dest="maxis",
	help="Create height (meters) level axis information"
)
parser.add_option(
	"--time-bounds", dest="tbounds", metavar="H1,H2",
	help="Create a time_bnds variable to specify the period of time considered in each time record. H1 is the start time in hours from the current time record and H2 is the ending time"
)
parser.add_option(
	"-r", "--reference-date", dest="refdate",
	help="Reference date for the files"
)
parser.add_option(
	"--time-units", dest="time_units", default="hours since 1950-01-01_00:00:00",
	help="Units for the time axis", metavar="Days/Hours since YYYY-MM-DD_hh:mm:ss (T or space as separator also work)"
)
parser.add_option(
	"-o", "--output", dest="OFILE", metavar="OUTPUTFILE.nc",
	help="Output file name"
)
parser.add_option(
	"--output-pattern", dest="OUTPUT_PATTERN", metavar="[varcf]_[varwrf]_[level]_[firsttime]_[lasttime]_experiment.nc",
	help="Output pattern to use if the option --split-in-variables is activated. Patterns recognized are currently of the form '[varcf]_[varwrf]_[firsttime]_[lasttime]_experiment.nc. Firsttime and lasttime are replaced by datetimes of the form YYYYmmddHH"
)
parser.add_option(
	"-g", "--geofile", metavar="geo_em.d0X.nc", dest="geofile",
	help="geo_em file to be used. For instance if you already removed geographic variables from the input files"
)
parser.add_option(
	"--fullfile", metavar="wrfout_d0X_allvars.nc", dest="fullfile",
	help="wrfout file to be used for variables not found in the data files. For instance if you removed variables as the eta or soil levels."
)
parser.add_option(
	"--split-variables", action="store_true", dest="splitvars",
	help="Write a separated file for each variable."
)
parser.add_option(
	"--split-levels", action="store_true", dest="splitlevs",
	help="Write a separated file for each variable and level. (Use always with --split-variables)"
)
parser.add_option(
	"--temp-dir", dest="TEMPDIR",
	help="Temporary directory to write the files while running. They are copied to the folfer specified by -o or --output-pattern when the program finishes."
)
parser.add_option(
	"--plevs-filter", dest="selected_plevs",
	help="Comma separated list of the pressure levels that we want to save, in hPa."
)
parser.add_option(
	"--slevs-filter", dest="selected_slevs",
	help="Comma separated list of the soil levels that we want to save, in meters. They must be in .2f format."
)
parser.add_option(
	"--output-format", dest="oformat", default="NETCDF3_CLASSIC",
	help="Format of the output files. Available possibilities: NETCDF4_CLASSIC, NETCDF3 (default). If using NETCDF4_CLASSIC, the deflate level is 4 by default."
)
parser.add_option(
	"--filter-times", dest="ftimes", default=False, metavar="%Yi%mi%di%Hi,%Yf%mf%df%Hf",
	help="Filter the output files so only times between the two selected are retained."
)
(opt, args) = parser.parse_args()
