#!/bin/bash

dir_id=$1

CHUNK_SIZE=100
NUM_FILES=$(expr `(ls -l gses_${dir_id}/matrix* | wc -l)` - 1)
NUM_JOBS=$((NUM_FILES / CHUNK_SIZE +1))

echo $NUM_JOBS
sbatch --array=1-${NUM_JOBS} submit_convert_to_obj.sh $dir_id 
