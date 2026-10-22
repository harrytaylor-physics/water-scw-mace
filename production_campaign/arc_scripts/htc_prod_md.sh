#!/bin/bash
#SBATCH --job-name=prod_md
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --array=0-53%6
#SBATCH --output=logs/prod_%A_%a.out
module load Anaconda3/2025.06-1
cd /home/shug7909/scw
line=$(tail -n +2 phase0_grid_v3.csv | sed -n "$((SLURM_ARRAY_TASK_ID+1))p")
T=$(echo "$line" | cut -d, -f1); RHO=$(echo "$line" | cut -d, -f2)
echo "task $SLURM_ARRAY_TASK_ID -> T=$T rho=$RHO"
python run_prod_md.py "$T" "$RHO"
