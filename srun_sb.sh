#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=users
#SBATCH --exclusive
#SBATCH -J simwave
#SBATCH --mail-user=long.qu@kaust.edu.sa
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00
#SBATCH --output=log/job.%J.out
#SBATCH --error=log/job.%J.err
#SBATCH --hint=nomultithread
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=40
##SBATCH --constraint=hwperf
#SBATCH -w cn512-22-l

#OpenMP settings:
export OMP_NUM_THREADS=40
export OMP_PROC_BIND=true
export OMP_PLACES=threads

#run the application:

echo $HOSTNAME
more /proc/cpuinfo

#srun ./bin/modeling  --verbose --n1  512 --n2  512 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301
#srun ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301
#srun ./bin/modeling  --verbose --n1 2048 --n2 2048 --n3 512  --iter 200 --dshot 1 --first 1301 --last 1301

#MODELING
#srun ./bin/modeling --verbose --n1 2048 --n2 2048  --n3 512 --iter 506  --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 24

MODE="MEM"
echo $MODE

srun ./bin/rtm --verbose --n1  512 --n2  512  --n3 512 --iter 48 --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 12
srun ./bin/rtm --verbose --n1 1024 --n2 1024  --n3 512 --iter 48 --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 12
srun ./bin/rtm --verbose --n1 2048 --n2 2048  --n3 512 --iter 24  --mode 2  --dshot 1 --first 1301 --last 1301  --nbsnap 12

MODE="I/O"
echo $MODE

srun ./bin/rtm --verbose --n1  512 --n2  512  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 4
srun ./bin/rtm --verbose --n1 1024 --n2 1024  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 4
srun ./bin/rtm --verbose --n1 2048 --n2 2048  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 4

srun ./bin/rtm --verbose --n1  512 --n2  512  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 8
srun ./bin/rtm --verbose --n1 1024 --n2 1024  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 8
srun ./bin/rtm --verbose --n1 2048 --n2 2048  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 8

srun ./bin/rtm --verbose --n1  512 --n2  512  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 12
srun ./bin/rtm --verbose --n1 1024 --n2 1024  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 12
srun ./bin/rtm --verbose --n1 2048 --n2 2048  --n3 512 --iter 48 --mode 1  --dshot 1 --first 1301 --last 1301  --nbsnap 12

