#!/bin/bash
#SBATCH --job-name=sex_lab
#SBATCH --output=sex_lab_%A_%a.out
#SBATCH --error=sex_lab_%A_%a.err
#SBATCH --time=1:00:00 
#SBATCH --mem=10000
#SBATCH --array=1
#SBATCH --partition=rbaltman

Rscript study_sex_label.R $SLURM_ARRAY_TASK_ID $1
