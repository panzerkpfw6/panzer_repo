#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=32
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export CFLAGS=-xHost
export CXXFLAGS=-xHost
###********** MODULES & COMPILING *********###
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5    #intel compiler: newer is better  #intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load intel-oneapi/2021.4.0
#module load intel/2022.1.0
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load cmake
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install
###*** PARAMETERS TO ITERATE A SWEEP ON ***###
# Values to iter for num_threads | tgs | size |
size=512
timesteps=1000 ###???
num_wf=21
t_dim=3
mode=2
###*********** RUNNING SIMWAVE ************###
tgs=4
th_z=$(expr $tgs / 2)
th_y=2
th_x=1
export FIRST_TOUCH=1
numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c