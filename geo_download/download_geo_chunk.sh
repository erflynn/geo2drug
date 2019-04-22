#!/bin/bash
#SBATCH --job-name=download
#SBATCH --output=download_geo_%A_%a.out
#SBATCH --error=download_geo_%A_%a.err
#SBATCH --time=2:00:00 
#SBATCH --mem=8000
#SBATCH --partition=rbaltman

ml python/3.6.1

export PYTHONPATH=${PYTHONPATH}:/scratch/users/erflynn/applications/python_packages/lib/python3.6/site-packages/
export PATH=${PATH}://scratch/users/erflynn/applications/python_packages//bin/

i=$SLURM_ARRAY_TASK_ID
i2=`printf "%02d" $i`
python3 downloadGEO.py tmp/gse${i2}