#!/bin/bash
#SBATCH --job-name=dip_lda
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --array=0-0
#SBATCH --output=/data/phys-plasma-dft/shug7909/lda_pipeline/logs/dip_lda_%A_%a.out
module purge
module load CP2K/2023.1-foss-2023a
export CP2K_DATA_DIR=${CP2K_DATA_DIR:-$EBROOTCP2K/data}
export OMP_NUM_THREADS=1
mapfile -t INPS < /data/phys-plasma-dft/shug7909/lda_pipeline/diplabels/dip_lda.list
inp="${INPS[$SLURM_ARRAY_TASK_ID]}"
if [ -z "$inp" ]; then echo "no input"; exit 0; fi
cd "$(dirname "$inp")" || exit 1
if [ -f dip.out ] && grep -q "SCF run converged" dip.out; then echo "skip"; exit 0; fi
mpirun cp2k.psmp -i dip.inp -o dip.out
rm -f *.wfn* *RESTART* 2>/dev/null
