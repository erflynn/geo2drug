#!/bin/bash
#SBATCH --job-name=combine_mat
#SBATCH --output=logs/combine_mat_%A_%a.out
#SBATCH --error=logs/combine_mat_%A_%a.err
#SBATCH --time=1:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1

ORGANISM=$1
OUT_DIR=data/03_silver_std/${ORGANISM}/02_keep_labels/
mkdir -p $OUT_DIR
Rscript code/02_process/02a_combine_sex_lab.R $ORGANISM $SLURM_ARRAY_TASK_ID
