#!/bin/bash
#SBATCH --job-name=opac_dil --partition=short --gres=gpu:1
#SBATCH --cpus-per-task=8 --mem=48G --time=02:00:00
#SBATCH --array=0-8%4 --output=logs/opac_dil_%A_%a.out
module load Anaconda3/2025.06-1
read T RHO N < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" cells_dilute.list)
d=$DATA/scw_prod/T${T}_rho${RHO}_05fs
python mace_opacity.py "$d" auto "$N" 0.5 water_dipole_v3.model 20 1
