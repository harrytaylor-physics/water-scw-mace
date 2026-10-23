#!/bin/bash
#SBATCH --job-name=diptrain
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/diptrain_%j.out
module load Anaconda3/2025.06-1
export PATH=$HOME/.local/bin:$PATH
mace_run_train \
  --name=water_dipole_v3 \
  --model=AtomicDipolesMACE --E0s='{1:0.0,8:0.0}' \
  --train_file=dipole_train_v3.xyz \
  --valid_fraction=0.1 \
  --dipole_key=REF_dipole \
  --loss=dipole \
  --error_table=DipoleRMSE \
  --r_max=6.0 --num_interactions=2 --num_channels=128 --max_L=2 \
  --max_num_epochs=200 --batch_size=10 \
  --default_dtype=float32 --device=cuda --seed=1
