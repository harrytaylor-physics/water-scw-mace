#!/bin/bash
#SBATCH --job-name=aimd_label
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --array=0-0
#SBATCH --output=logs/aimd_label_%A_%a.out
module purge
module load CP2K/2023.1-foss-2023a
export CP2K_DATA_DIR=${CP2K_DATA_DIR:-$EBROOTCP2K/data}
export OMP_NUM_THREADS=1
mapfile -t INPS < frames_aimd.list
inp="${INPS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$inp" ]; then echo "no input for task $SLURM_ARRAY_TASK_ID"; exit 0; fi
dir=$(dirname "$inp"); cd "$dir" || exit 1
echo "[$SLURM_ARRAY_TASK_ID] $(hostname) -> $dir"
mpirun cp2k.psmp -i input.inp -o output.out
