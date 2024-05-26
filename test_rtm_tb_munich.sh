#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=test_rtm_2nd
#SBATCH --output=logs/rtm.%J.out
#SBATCH --error=logs/rtm.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-rtm.log) 2>&1
echo $hostname
#lscpu
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

###********** MODULES & COMPILING *********###
module load icc/2020.2.254
module load cmake

#####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
#timesteps=2000
timesteps=2000
nb_snap=100
export size1=2048
export size2=512
#export size=256
#####*********** SB tests ************###
echo !!SB!!
#numactl --interleave=all ./bin/rtm --verbose --n1 $size  --n2 $size --n3 $size --iter $timesteps --dshot 1 --first 1301 --last 1301
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --unbuffered numactl --interleave=all ./bin/rtm --verbose --n1 $size1  --n2 $size1 --n3 $size2 --iter $timesteps --dshot 1 --first 1301 --last 1301 --nbsnap $nb_snap
#####*********** SB tests ************###
echo !!TB default!!
#numactl --interleave=all ./bin/rtm --cpu --verbose --n1 $size  --n2 $size --n3 $size --iter $timesteps --dshot 1 --first 1301 --last 1301
#####*********** TB tests ************########***********
echo !!TB with parameters!!
#export FIRST_TOUCH=1m
#srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
mode=2
th_x=8
th_y=2
th_z=1
tgs=16
t_dim=7
num_wf=64
srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --unbuffered numactl --interleave=all ./bin/rtm --verbose --n1 $size1  --n2 $size1 --n3 $size2 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c --nbsnap $nb_snap

####*********** RUNNING SIMWAVE ************###
## iterate on parameters values
#echo "CSV: size ; timesteps ; OMP_NUM_THREADS ; tgs ; ntg ; (th_x, th_y, th_z) ; t_dim ; num_wf ; mode ; tb or sb"
#size=512
## Thread configuration: SIMWAVE x=1, RACHED z=1 (innermost dimension is vectorized)
####*** temporal blocking parameters***###
##thread group : (1,2,2)  #group size: 4 #num_group : 8
#th_z=2
#th_y=2
#th_x=1
#tgs=4
#num_wf=21
#t_dim=3
#echo $(expr $OMP_NUM_THREADS / $tgs)
#echo "PROGRAM CALL: numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c"
#echo "CSV: $size ; $timesteps ; $OMP_NUM_THREADS ; $tgs ; $(expr $OMP_NUM_THREADS / $tgs) ; ($th_x, $th_y, $th_z) ; $t_dim ; $num_wf ; $mode ; tb"
####*********** TB tests ************###
## Run without First Touch
#export FIRST_TOUCH=0
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
#
### Run with First Touch
#export FIRST_TOUCH=1
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
#####*********** SB tests ************###
##numactl --interleave=all ./bin/modeling-snap_end --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
#
#echo modelling_without_snapshot_saving
#echo SB
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
#echo SB_different
#srun --unbuffered numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##sleep 600
#
