#!/bin/bash
#SBATCH --job-name=train_v4
#SBATCH --partition=short
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=logs/train_v4_%j.out
module load Anaconda3/2025.06-1
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
cd /home/shug7909/scw
mace_run_train --restart_latest \
  --name=water_scw_v4 \
  --foundation_model=/home/shug7909/scw/mace-mpa-0-medium.model \
  --multiheads_finetuning=False \
  --train_file=train_v4.xyz \
  --energy_key=energy --forces_key=forces \
  --E0s=average \
  --energy_weight=1 --forces_weight=100 \
  --loss=weighted \
  --swa --start_swa=120 \
  --max_num_epochs=200 \
  --batch_size=2 --valid_batch_size=2 \
  --default_dtype=float32 --device=cuda \
  --seed=1
