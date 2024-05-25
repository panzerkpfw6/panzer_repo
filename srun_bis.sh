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


#OpenMP settings:
export OMP_NUM_THREADS=40
export OMP_PROC_BIND=true
export OMP_PLACES=cores

#run the application:

echo $HOSTNAME
more /proc/meminfo
srun numactl --interleave=all ./bin/modeling  --verbose --n1 2048 --n2 2048 --n3 512  --iter 88  --tb_thread_group_size 10 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15   --mode 3  --src_depth 50 --rcv_depth 50  --dshot 1 --first 0 --last 0  --fwd_steps 3 -c

srun numactl --interleave=all ./bin/rtm  --verbose --n1 2048 --n2 2048 --n3 512  --iter 88  --tb_thread_group_size 10 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15   --mode 2  --dshot 1 --first  0 --last 0  --fwd_steps 3 -c


