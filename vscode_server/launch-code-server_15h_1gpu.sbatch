#!/bin/bash --login
#SBATCH --time=13:55:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-gpu=6
#SBATCH --mem=64G
#SBATCH --partition=batch
#SBATCH --reservation=A100
#SBATCH --job-name=launch-code-server
#SBATCH --mail-type=ALL
#SBATCH --output=./%x-%j-slurm.out
#SBATCH --error=./%x-%j-slurm.err

# use srun to launch code server in order to reserve a port
srun --resv-ports=1 ./launch-code-server.srun
