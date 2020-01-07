#!/bin/bash
#SBATCH --job-name=exp_mat
#SBATCH --output=exp_mat_%A_%a.out
#SBATCH --error=exp_mat_%A_%a.err
#SBATCH --time=4:00:00 
#SBATCH --mem=10000
#SBATCH --array=2-12
#SBATCH --partition=rbaltman

Rscript create_exp_mat.R $SLURM_ARRAY_TASK_ID
