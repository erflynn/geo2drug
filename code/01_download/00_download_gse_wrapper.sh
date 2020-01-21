#!/bin/bash
# download_gse_wrapper.sh
# E Flynn
# 4/22/2019
#
# wrapper for downloading a lot of GSEs in parallel
# takes as input a list of GSEs, downloads all to a specified directory.

ID=$1
gse_list=$2  #'gse_for_silver_std.csv'

mkdir -p logs
mkdir -p tmp_${ID}
mkdir -p gses_${ID}
mkdir -p gses_${ID}/matrix
split -l 100 -d ${gse_list} tmp_${ID}/gse
NUM_JOBS=$(expr `(ls -l tmp_${ID}/gse* | wc -l)` - 1)
echo $NUM_JOBS
sbatch --array=0-${NUM_JOBS} code/01_download/00a_download_geo_chunk.sh ${ID}


# after
# grep "ERROR" *.out > error_gses.txt
