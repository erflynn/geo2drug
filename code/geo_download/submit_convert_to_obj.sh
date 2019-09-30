#!/bin/bash
#SBATCH --job-name=check_geo
#SBATCH --error=logs/check_geo_%A_%a.out
#SBATCH --error=logs/check_geo_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

dir_id=$1

mkdir -p gses_${dir_id}/rObj
Rscript seriesToObj.R $SLURM_ARRAY_TASK_ID $dir_id
