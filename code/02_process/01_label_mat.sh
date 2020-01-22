#!/bin/bash
#SBATCH --job-name=label_mat
#SBATCH --output=logs/label_mat_%A_%a.out
#SBATCH --error=logs/label_mat_%A_%a.err
#SBATCH --time=5:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1

ORGANISM=$1
OUT_DIR=data/03_silver_std/${ORGANISM}/01_compare_labels/
mkdir -p $OUT_DIR
Rscript code/02_process/01a_study_sex_label.R $OUT_DIR $ORGANISM $SLURM_ARRAY_TASK_ID
