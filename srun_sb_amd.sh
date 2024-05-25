#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=log/job.%J.out
#SBATCH --error=log/job.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=long.qu@kaust.edu.sa
#SBATCH --time=24:00:00
#SBATCH --partition=users
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH -w cn110-22-l

echo $HOSTNAME

#OpenMP settings:
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=true
export OMP_PLACES=cores

#run the application:

echo $HOSTNAME
more /proc/cpuinfo

#srun ./bin/modeling  --verbose --n1  512 --n2  512 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301
#srun ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301
#srun ./bin/modeling  --verbose --n1 2048 --n2 2048 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301

srun ./bin/rtm --verbose --n1  512 --n2  512  --n3 512 --iter 96 --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 24
srun ./bin/rtm --verbose --n1 1024 --n2 1024  --n3 512 --iter 96 --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 24
srun ./bin/rtm --verbose --n1 2048 --n2 2048  --n3 512 --iter 96 --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 24

