#!/bin/bash
#SBATCH --job-name=prod_ext --partition=medium --gres=gpu:1
#SBATCH --cpus-per-task=8 --mem=64G --time=24:00:00
#SBATCH --array=0-20%6 --output=logs/prod_ext_%A_%a.out
module load Anaconda3/2025.06-1
read T RHO N < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cells_ext.list)
python run_prod_md_ext.py "$T" "$RHO" "$N"
