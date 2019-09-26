#!/bin/bash
#SBATCH --job-name=check_geo
#SBATCH --output=check_geo_%A_%a.out
#SBATCH --error=check_geo_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=8000
#SBATCH --array=1,2,4-12
#SBATCH --partition=rbaltman

Rscript seriesToObj.R $SLURM_ARRAY_TASK_ID
