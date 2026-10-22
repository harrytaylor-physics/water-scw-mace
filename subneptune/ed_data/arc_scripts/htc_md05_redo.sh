#!/bin/bash
#SBATCH --job-name=md05r --partition=short --gres=gpu:1
#SBATCH --cpus-per-task=8 --mem=32G --time=03:00:00
#SBATCH --output=logs/md05redo_%j.out
module load Anaconda3/2025.06-1
cd /home/shug7909/scw
python run_prod_md_05fs.py 3000 108.796
