#!/bin/bash
#SBATCH --job-name=opac05 --partition=short --gres=gpu:1
#SBATCH --cpus-per-task=8 --mem=32G --time=01:00:00 --array=0-14%4
#SBATCH --output=logs/opac05_%A_%a.out
module load Anaconda3/2025.06-1
cd /home/shug7909/scw
mapfile -t D < cells_05_dirs.list
d="${D[$SLURM_ARRAY_TASK_ID]}"
[ -z "$d" ] && exit 0
python mace_opacity.py "$d" auto 64 0.5 water_dipole_v3.model 20 1
