#!/bin/bash
#SBATCH --job-name=train_run
#SBATCH --output=logs/train_run_%A.out
#SBATCH --error=logs/train_run_%A.err
#SBATCH --time=2:00:00 
#SBATCH --mem=10000
#SBATCH --partition=rbaltman

ml R/3.6.1

ORGANISM=$1
RUN_TYPE=$2
DAT_FILE=$3 
# for test_dat: "data/01_sample_lists/<org>_testing_<run_type>.csv", etc (no suffix for common)
# for train_dat: "data/01_sample_lists/<org>_training_<run_type>.csv"
SUB_DIR=$4 # `03_silver_std/` or `05_single_sex/`

MY_DIR=data/${SUB_DIR}/${ORGANISM}
mkdir -p ${MY_DIR}/03_out_mat/
Rscript code/02_process/05d_construct_mat.R $ORGANISM $SLURM_ARRAY_TASK_ID $RUN_TYPE $DAT_FILE $MY_DIR