#!/bin/bash

function download(){
  url=$1
  outfile=$(basename ${url})
  test -f ${outfile} || wget ${url}
}

# Get test data from the WRF tutorial
download https://www2.mmm.ucar.edu/wrf/TUTORIAL_DATA/single_domain/geo_em.d01.nc
download https://www2.mmm.ucar.edu/wrf/TUTORIAL_DATA/single_domain/wrfout_d01.tar.gz
tar xzvf wrfout_d01.tar.gz

# Extract orography from geo_em
wrfncxnj -g geo_em.d01.nc --single-rec -v HGT_M -o test_orog.nc 

# Process pr,tas and huss
wrfncxnj -v RAINF,T2,Q2 --split-variables --output-pattern='test_[varcf]_[firsttime]_[lasttime].nc' wrfout_d01_2016-10-*

# Process 3D data
# Requires custom table with 'z' level type instead of the default 'p', which assumes pressure level interpolation through p_interp
#wrfncxnj -v QVAPOR -t wrfncxnj.table --split-variables --output-pattern='test_[varcf]_[firsttime]_[lasttime].nc' wrfout_d01_2016-10-*

# Process soil data
wrfncxnj -v SMOIS --split-variables --output-pattern='test_[varcf]_[firsttime]_[lasttime].nc' wrfout_d01_2016-10-*
