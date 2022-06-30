
## wrfncxnj - WRF NetCdf eXtract aNd Join


The output of the Weather Research and Forecasting model (WRF) is, by default, 
written in NetCDF files formatted according to specific WRF conventions. Most 
of the software packages (e.g. CDO) and international projects that involve 
climate modelling (CORDEX, CMIP, etc.) require CF conventions. Also, the way 
output data are produced by the model in a climate run (all variables in a 
single file split in time) does not match the modeller's usage pattern (usually
requiring a few variables over the whole period). WRF NetCDF Extract&Join 
(wrfncxnj) addresses this problem in a flexible way. It is written in python, 
but is its designed to be invoked from the command line. Due to the large 
amount of available options it is usually preferable to call wrfncxnj from a 
Shell script. The main goal is to simplify the post-processing of large amounts
of data generated by WRF, making it a less time-consuming and error-prone 
task. Since python is an interpreted language, WRF NetCDF Extract&Join does 
not need to be installed or compiled, and it has few dependencies, that can 
be installed from repositories in most linux distributions.

The tool has been originally developed at the University of Cantabria by 
Markel García-Díez, Jesús Fernández, Lluís Fita and Álvaro Lavin


### Installation

```
python setup.py install
```


It has the following system dependencies:

* NetCDF libraries with NetCDF4 support http://www.unidata.ucar.edu/software/netcdf.

* Python http://www.python.org/

* Numpy http://numpy.scipy.org/

* NetCDF4 python library http://unidata.github.io/netcdf4-python/

Note that the netCDF4 python library is not equivalent to the python-netcdf 
package available from linux repositories. Currently, it must be downloaded 
from the web cited above.

Another possibility (easier) is to install a scientific python 
distribution such as [Anaconda https://store.continuum.io/cshop/anaconda/], 
that already comes with this packages. 

### Usage

WRF NetCDF Extract&Join is invoked from the command line as:

```
wrfncxnj [options] [files to process]
```

A description of all the options available can be obtained executing:

```
wrfncxnj.py -h
```

The ASCII file wrfncxnj.table plays a central role in translating the WRF 
variables to CF conventions. This file can be modified by the user to meet 
its needs. The user interested in modifying the rest of the code should check 
the /TechnicalDescription. Users are highly encouraged to look into the code, 
play with it, and report back any bug or improvement.

The essential elements of the command line needed are the list of variables 
to be extracted, the output file (or pattern, is the output is to be splitted 
into variables or levels), the reference time for the time axis and the 
"wrfout" files, which can be provided as arguments. See the following example:

```
wrfncxnj.py [flags] --split-variables --time-units="days since 1949-12-01_00:00:00" --output-pattern=EXPERIMENT_[varcf].nc -v "VAR1,VAR2,VAR3" wrfout1.nc wrfout2.nc ... wrfoutN.nc
```

This command would read variables VAR1, VAR2 and VAR3 from wrfout1.nc 
wrfout2.nc ... wrfoutN.nc and write them in separate files labelled by their 
CF name. Notice that VAR1, VAR2 and VAR3 are not the CF names of the 
variables, but the name that appears in the first column of wrfncxnj.table 
(usually, the variable name in WRF). Different wrfncxnj variables can share 
the same CF name, so the unique names for the variables are the WRF names 
(or the aliases in the first column of wrfncxnj.table).