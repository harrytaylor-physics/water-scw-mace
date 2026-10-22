#!/bin/bash
#SBATCH --job-name=prod_dil --partition=medium --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8 --mem=64G --time=24:00:00
#SBATCH --array=0-8%4 --output=logs/prod_dil_%A_%a.out
module load Anaconda3/2025.06-1
read T RHO N < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cells_dilute.list)
python run_prod_md_ext.py "$T" "$RHO" "$N"
