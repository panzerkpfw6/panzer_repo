#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=workq
#SBATCH --job-name=test_default_pars
#SBATCH --output=logs/rtm_salt3d_tb.%J.out
#SBATCH --error=logs/rtm_salt3d_tb.%J.err
#SBATCH --cpus-per-task=192
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** COMMENT *********###
# In this test we check performance with best SB,TB parameters for AMD MilanX.
# be aware that cache blocking should be different for 1st order and second order
# 1)please change SLURM CONFIGURATION accordingly
#  partition,
# number of cpus-per-task should be equal to number of available cores.
# 2)please change OPENMP PARAMETERS accordingly
#OMP_NUM_THREADS should be equal to number of cpus
# 3)please change MODULES accordingly
#we need to load cmake,icpc modules from somewhere

###******** HORODATED LOG WRITING *********###
mkdir velocity_models
cd velocity_models
#wget -q --content-disposition 'https://www.dropbox.com/scl/fi/05bc42ctyyhx1d7d10k51/salt3d_676x676x201_xyz.raw?rlkey=04646f0r2vil9ph5m9yjsiras&dl=0'
cd ..

echo $hostname
lscpu

rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
#rm ./data/*ilm*
#rm ./data/*img*
#rm ./data/*sismos*
#rm ./data/*snap*

###********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=192
##export OMP_NUM_THREADS=64
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
#export OMP_NESTED='True'
#export granularity=fine
#export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t

###********** Set compiler flags *********###
#export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

######################################################
export OMP_PLACES=cores;
export OMP_PROC_BIND=close;
export OMP_STACKSIZE=64M;
export OMP_NUM_THREADS=192;
export CFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export CXXFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export FFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"

###********** MODULES *********###
########module load intel/2024.2.1
module load intel-oneapi/2023.1.0
module load cmake

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
rm -rf ./logs/rtm_salt_tb #delete if existing
mkdir ./logs/rtm_salt_tb
export logs_path=./logs/rtm_salt_tb

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
timesteps=4001; dt=0.001;
nx=672;ny=672;nz=201;dh=25;
##### Profile x=310. Salt3D. no mistake
#first=20957;last=21023;
#### Profile x=?. Salt3D. no mistake
first=1;last=98;
#first=1;last=3;
dshot=4568;
fmax=11;
cbx=8;cby=6;cbz=9999;
##################  TB parameters  ###################
x=3; y=2; z=2; t=3; w=24;
tgs=$(expr $x \* $y \* $z); echo $tgs
######################################################

#echo "Model data for RTM. salt3d."
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
#--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
#--tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#--tb_t_dim $t --tb_num_wf $w --dx $dh --dy $dh --dz $dh --dt $dt -c >> $logs_path/log_model_salt3d.log;
#
#echo "Perform RTM"
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
#--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
#--tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#--tb_t_dim $t --tb_num_wf $w --dx $dh --dy $dh --dz $dh --dt $dt -c >> $logs_path/log_rtm_salt3d.log;

echo "Gather images"
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
--first $first --last $last -c --fwd_steps 3 --order 1 --src_depth 5 --rcv_depth 8 --drcv 1 --dir "./data" \
--tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
--tb_t_dim $t --tb_num_wf $w --dx $dh --dy $dh --dz $dh --dt $dt >> $logs_path/log_gather_salt3d.log;


