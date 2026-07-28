#!/bin/bash
#SBATCH --job-name=lda_anchor
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --account=phys-plasma-dft
#SBATCH -o /data/phys-plasma-dft/shug7909/lda_pipeline/slurm-%j.out
module load CP2K/2023.1-foss-2023a
export OMP_NUM_THREADS=1
for INP in /data/phys-plasma-dft/shug7909/lda_pipeline/anchor/T*_rho*/anchor_*/label.inp; do
  DIR=$(dirname "$INP")
  if [ -f "$DIR/label.out" ] && grep -q "SCF run converged" "$DIR/label.out"; then
    echo "skip: $DIR"; continue
  fi
  cd "$DIR"
  mpirun -np 48 --bind-to core cp2k.psmp -i label.inp -o label.out
done
