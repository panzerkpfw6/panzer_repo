#!/bin/bash
#SBATCH -N 1
#SBATCH --partition=users
#SBATCH --exclusive
#SBATCH -J stencil
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
#export OMP_PROC_BIND=true
#export OMP_PLACES=cores

#run the application:

echo $HOSTNAME



echo $TESTNAME

export KMP_AFFINITY=verbose,granularity=core,scatter

export OMP_NUM_THREADS=40
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 

export OMP_NUM_THREADS=30
export KMP_AFFINITY="granularity=fine,verbose,proclist=[0,1,2,4,5,6,8,9,10,12,13,14,16,17,18,20,21,2
2,24,25,26,28,29,30,32,33,34,36,37,38],explicit"
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 3 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 
export OMP_NUM_THREADS=20
export KMP_AFFINITY="granularity=fine,verbose,proclist=[0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32
,33,36,37],explicit"
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 2 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 --tb_affinity affinity2_csl.txt

export OMP_NUM_THREADS=10
export KMP_AFFINITY="granularity=fine,verbose,proclist=[0,4,8,12,16,20,24,28,32,36],explicit"
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 10 --tb_nb_thread_groups 1 --tb_th_x 1 --tb_th_y 2 --tb_th_z 5 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 --tb_affinity affinity2_csl.txt

export OMP_NUM_THREADS=1
srun numactl --interleave=all ./bin/modeling  --verbose --n1 1024 --n2 1024 --n3 512  --iter 504  --tb_thread_group_size 1 --tb_nb_thread_groups 1 --tb_th_x 1 --tb_th_y 1 --tb_th_z 1 --tb_t_dim 3 --tb_num_wf 15    --mode 2  --src_depth 50 --rcv_depth 50  --dshot 1 --first 1301 --last 1301  --fwd_steps 4 
