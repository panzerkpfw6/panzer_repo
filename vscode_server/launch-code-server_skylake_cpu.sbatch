#!/bin/bash --login
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --partition=batch
#SBATCH --constraint=[skylake]
#SBATCH --mem=32G
####SBATCH --reservation=A100

#SBATCH --job-name=launch-cs-skylake
#SBATCH --mail-type=ALL
#SBATCH --output=./%x-%j-slurm.out
#SBATCH --error=./%x-%j-slurm.err

# use srun to launch code server in order to reserve a port
srun --resv-ports=1 ./launch-code-server.srun
