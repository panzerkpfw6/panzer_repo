#!/bin/bash --login
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
##SBATCH --cpus-per-gpu=6  
##SBATCH --constraint=[v100]
##SBATCH --constraint=[rtx2080ti]
##SBATCH --constraint=[gtx1080ti]
##SBATCH --constraint=[p100]
##SBATCH --gres=gpu:rtx2080ti:1
#SBATCH --mem=16G
#SBATCH --partition=debug 
##SBATCH --partition=batch
#SBATCH --job-name=launch-code-server
#SBATCH --mail-type=ALL
#SBATCH --output=./%x-%j-slurm.out
#SBATCH --error=./%x-%j-slurm.err

# use srun to launch code server in order to reserve a port
# source ~/.bashrc_ibex
srun --resv-ports=1 ./launch-code-server.srun
