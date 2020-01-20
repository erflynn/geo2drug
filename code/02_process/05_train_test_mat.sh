#!/bin/bash
#SBATCH --job-name=train_run
#SBATCH --output=logs/train_run_%A.out
#SBATCH --error=logs/train_run_%A.err
#SBATCH --time=2:00:00 
#SBATCH --mem=20000
#SBATCH --partition=rbaltman

ml R/3.6.1

Rscript code/02_process/05a_train_test_mat.R
