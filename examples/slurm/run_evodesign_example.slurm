#!/bin/bash
#SBATCH --job-name=evodesign
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --output=evodesign.out
#SBATCH --error=evodesign.err

cd $SLURM_SUBMIT_DIR

module load python/miniforge
source ~/.bashrc
conda activate evodesign
./run_from_terminal.sh
