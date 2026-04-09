#!/bin/bash --login
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --partition=rome
#SBATCH --mem=32G

#SBATCH --job-name=launch-code-server
#SBATCH --mail-type=ALL
#SBATCH --output=./%x-%j-slurm.out
#SBATCH --error=./%x-%j-slurm.err

# use srun to launch code server in order to reserve a port
srun --resv-ports=1 ./launch-code-server_ibex_cpu.srun
