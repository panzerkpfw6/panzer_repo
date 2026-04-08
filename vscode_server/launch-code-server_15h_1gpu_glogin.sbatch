#!/bin/bash --login
#SBATCH --time=15:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=6  
##SBATCH --reservation=v10
#SBATCH --gres=gpu:1
#SBATCH --constraint=[v100]
#SBATCH --mem=64G
##SBATCH --partition=debug 
#SBATCH --partition=batch
#SBATCH --job-name=launch-code-server
#SBATCH --mail-type=ALL
#SBATCH --output=./%x-%j-slurm.out
# use srun to launch code server in order to reserve a port
# source ~/.bashrc_ibex
srun --resv-ports=1 ./launch-code-server.srun
