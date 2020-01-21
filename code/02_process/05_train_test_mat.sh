#!/bin/bash
#SBATCH --job-name=train_run
#SBATCH --output=logs/train_run_%A.out
#SBATCH --error=logs/train_run_%A.err
#SBATCH --time=2:00:00 
#SBATCH --mem=10000
#SBATCH --partition=rbaltman

ml R/3.6.1

ORGANISM="human"
Rscript code/02_process/05b_test_mat.R ${ORGANISM} $SLURM_ARRAY_TASK_ID
