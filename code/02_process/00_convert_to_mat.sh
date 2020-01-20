#!/bin/bash
#SBATCH --job-name=convert_mat
#SBATCH --output=logs/convert_mat_%A_%a.out
#SBATCH --error=logs/convert_mat_%A_%a.err
#SBATCH --time=5:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1

DATA_DIR="/scratch/users/erflynn/sex_labeling/geo_pipeline/data"
GSE_FILE=${DATA_DIR}/sample_lists/gse_for_silver_std_human.csv
OUT_DIR=${DATA_DIR}/silver_std/
Rscript gse_to_mat.R $GSE_FILE $OUT_DIR $SLURM_ARRAY_TASK_ID
