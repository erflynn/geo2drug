#!/bin/bash
#SBATCH --job-name=convert_mat
#SBATCH --output=logs/convert_mat_%A_%a.out
#SBATCH --error=logs/convert_mat_%A_%a.err
#SBATCH --time=5:00:00 
#SBATCH --mem=15000
#SBATCH --partition=rbaltman

ml R/3.6.1
ORGANISM=$1
DS_DIR=$2
GSE_FILE=$3
GSE_PATH=data/01_sample_lists/${GSE_FILE}

OUT_DIR=data/${DS_DIR}/${ORGANISM}/00_mat_files/
mkdir -p ${OUT_DIR}
Rscript code/02_process/00a_gse_to_mat.R $GSE_PATH $OUT_DIR $ORGANISM $SLURM_ARRAY_TASK_ID

