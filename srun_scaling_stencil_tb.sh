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
##SBATCH --hint=nomultithread
##SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=40
##SBATCH -w cn512-22-l
#SBATCH --constraint=[cascadelake]
#SBATCH --reservation=testing-slurm-20.11

#OpenMP settings:
export OMP_NUM_THREADS=40
export OMP_PROC_BIND=true
export OMP_PLACES=cores

#run the application:

echo $HOSTNAME



echo $TESTNAME

#srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 604  --tb_thread_group_size 10 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c

xsrun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 3 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 2 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c --tb_affinity affinity2_csl.txt
srun numactl --interleave=all ./bin/modeling --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 1 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c --tb_affinity affinity2_csl.txt
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 1 --tb_nb_thread_groups 1 --tb_th_x 1 --tb_th_y 1 --tb_th_z 1 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 -c






