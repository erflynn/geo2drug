#!/bin/bash
#SBATCH --job-name=combine_mat
#SBATCH --output=logs/combine_mat_%A_%a.out
#SBATCH --error=logs/combine_mat_%A_%a.err
#SBATCH --time=1:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml R/3.6.1

DATA_DIR="/scratch/users/erflynn/sex_labeling/geo_pipeline/data"
IN_DIR=${DATA_DIR}/compare_labels/
cd $IN_DIR;
Rscript /scratch/users/erflynn/sex_labeling/geo_pipeline/code/combine_sex_lab.R $SLURM_ARRAY_TASK_ID
