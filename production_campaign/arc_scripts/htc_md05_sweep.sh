#!/bin/bash
#SBATCH --job-name=md05 --partition=short --gres=gpu:1
#SBATCH --cpus-per-task=8 --mem=32G --time=04:00:00 --array=0-14%4
#SBATCH --output=logs/md05_%A_%a.out
module load Anaconda3/2025.06-1
cd /home/shug7909/scw
mapfile -t PAIRS < cells_05.list
read T rho <<< "${PAIRS[$SLURM_ARRAY_TASK_ID]}"
echo "cell T=$T rho=$rho"
python run_prod_md_05fs.py "$T" "$rho"
