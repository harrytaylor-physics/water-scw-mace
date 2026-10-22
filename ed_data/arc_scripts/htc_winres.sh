#!/bin/bash -l
#SBATCH --job-name=winres --partition=short
#SBATCH --cpus-per-task=8 --mem=16G --time=00:30:00
#SBATCH --export=NONE
#SBATCH --output=logs/winres_%j.out
module purge
module load Anaconda3/2025.06-1
unset PYTHONPATH PYTHONNOUSERSITE
export OMP_NUM_THREADS=8
cd ~/scw && python run_winres_cpu.py
