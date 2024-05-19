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
module purge
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load cuda
module load cmake
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
#make clean
#make VERBOSE=1
#make install
###*** PARAMETERS TO ITERATE A SWEEP ON ***###
# Values to iter for num_threads | tgs | size |
size=256
num_wf=20
t_dim=7
mode=2
###*********** RUNNING SIMWAVE ************###
tgs=4
th_x=2
th_y=2
th_z=1
export FIRST_TOUCH=1
echo $OMP_NUM_THREADS
echo $tgs
echo $th_x
echo $th_y
echo $th_z
echo $t_dim
echo $num_wf
echo $mode
./bin/modeling --verbose --n1 256 --n2 256 --n3 256 --iter 505 --tb_thread_group_size 4 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
#numactl --interleave=all ./bin/modeling --verbose --n1 256 --n2 256 --n3 256 --iter 505 --tb_thread_group_size 4 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c