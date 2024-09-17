#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --nodelist=gbt04
#SBATCH --time=72:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=rtm
#SBATCH --output=./logs/rtm.%J.out
#SBATCH --error=./logs/rtm.%J.err
#SBATCH --hint=nomultithread

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-rtm.log) 2>&1
echo $hostname
#lscpu

####**********  Munich ***********###
#OpenMP settings:
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t
module load icc/2020.2.254
module load cmake

mkdir ./logs
mkdir ./logs/sb_param_search
export logs_path=./logs/sb_param_search/$SLURM_JOB_PARTITION
mkdir $logs_path
exec > >(while read line; do echo "$(date): $line"; done | tee ./logs/rtm_perf_$SLURM_JOB_PARTITION.log) 2>&1

echo "compilation"
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
#timesteps=2000
timesteps=2200
nt=2200;nt2=2*nt;
nx=128;ny=256;nz=512;
nx=128;ny=256;nz=128;
nx=1024;ny=1024;nz=512;
first=1;last=10;
first=16449;last=16449;
dshot=300;fmax=11;
###********** TB parameters ***********###
#num_th=128;th_x=8;th_y=2;th_z=1;num_wf=64;tdim=7; #stencil parameters. rewrite to simwave convention.
num_th=128;th_x=1;th_y=2;th_z=8;num_wf=64;tdim=7;tgs=th_x*th_y*th_z;
echo "${logs_path}/log-TB_1st_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}.log"

####*********** Modeling ************###
srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1 >> "${logs_path}/log_mod_SB_1st.log"
srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt2 --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $tdim --tb_num_wf $num_wf -c --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1 >> "${logs_path}/log_mod_TB_1st.log"
####*********** RTM ************###
srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1 >> "${logs_path}/prepare_sismos_rtm.log"
srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1 >> "${logs_path}/log_rtm_SB_1st.log"
srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $tdim --tb_num_wf $num_wf -c --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1 >> "${logs_path}/log_rtm_TB_1st.log"
