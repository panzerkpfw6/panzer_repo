#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=48
#export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t
###********** MODULES *********###
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load icc/2020.2.254
module load cmake

##### COMPILATION #####
#mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
####CC=icc CXX=icpc -DCMAKE_BUILD_TYPE=Debug cmake .
#CC=icc CXX=icpc cmake .
#make clean
#make VERBOSE=1
#make install

###############################
nx=128;ny=256;nz=512;
nt=57; dt=0.001;
shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
src_depth=256;
x=2; y=2; z=1; t=7; w=20; tgs=4;
export OMP_NUM_THREADS=4

./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
 --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot  --fwd_steps 3 -c --src_depth 256 --order 1 --fmax 8;
./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
 --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot  --fwd_steps 3 -c --src_depth 256 --order 2 --fmax 8;