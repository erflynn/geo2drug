#!/bin/bash
#SBATCH --job-name=label_mat
#SBATCH --output=logs/label_mat_%A_%a.out
#SBATCH --error=logs/label_mat_%A_%a.err
#SBATCH --time=5:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1


GSE_FILE=data/01_sample_lists/gse_for_silver_std_human.csv
OUT_DIR=data/03_silver_std/02_compare_labels/
Rscript code/02_process/01a_study_sex_label.R $GSE_FILE $OUT_DIR  $SLURM_ARRAY_TASK_ID
