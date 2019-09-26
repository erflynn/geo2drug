#!/bin/bash
# download_gse_wrapper.sh
# E Flynn
# 4/22/2019
#
# wrapper for downloading a lot of GSEs in parallel
# takes as input a list of GSEs, downloads all to a specified directory

mkdir -p tmp
mkdir -p gses
mkdir -p gses/matrix
gse_list=$1  #'gse_for_silver_std_gse.csv'
split -l 100 -d ${gse_list} tmp/gse
NUM_JOBS=$(expr `(ls -l tmp/gse* | wc -l)` - 1)
sbatch --array=0-${NUM_JOBS} download_geo_chunk.sh


# after
#rm -rf tmp
#grep "ERROR" *.out > error_gses.txt
#rm *.out
#rm *.err