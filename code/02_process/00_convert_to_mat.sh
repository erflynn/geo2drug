#!/bin/bash
#SBATCH --job-name=convert_mat
#SBATCH --output=logs/convert_mat_%A_%a.out
#SBATCH --error=logs/convert_mat_%A_%a.err
#SBATCH --time=5:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1
ORGANISM=$1
GSE_FILE=data/01_sample_lists/gse_for_silver_std_${ORGANISM}.csv
if [ ! -e "$GSE_FILE" ]; then 
    echo "Please specify another organism"
else
    OUT_DIR=data/03_silver_std/${ORGANISM}/00_mat_files/
    mkdir -p ${OUT_DIR}
    Rscript code/02_process/00a_gse_to_mat.R $GSE_FILE $OUT_DIR $ORGANISM $SLURM_ARRAY_TASK_ID
fi
