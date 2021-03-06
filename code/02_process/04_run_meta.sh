#!/bin/bash
#SBATCH --job-name=meta_run
#SBATCH --output=logs/meta_run_%A.out
#SBATCH --error=logs/meta_run_%A.err
#SBATCH --time=5:00:00 
#SBATCH --mem=20000
#SBATCH --partition=rbaltman

ml R/3.6.1

ORGANISM=$1
RUN_TYPE=$2 # common or full
mkdir -p data/03_silver_std/${ORGANISM}/04_meta_res
Rscript code/02_process/04a_run_meta.R $ORGANISM $RUN_TYPE
