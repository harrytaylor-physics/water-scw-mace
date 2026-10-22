#!/bin/bash
#SBATCH --job-name=spec05
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/spec05_%j.out
module load Anaconda3/2025.06-1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd /home/shug7909/scw
python run_bulk_md_05fs.py && python dipole_eval_05.py && python spectrum_05.py
