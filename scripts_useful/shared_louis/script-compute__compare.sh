#!/bin/bash
###********* SLURM CONFIGURATION **********###
#

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-compute__compare.log) 2>&1

###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=56
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'

###********** MODULES & COMPILING *********###
#intel compiler needed: newer is better
#last CMakeCache.txt is saved
module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
module load cmake
mv -f ./CMakeCache.txt ./CMakeCache-old.txt
export CFLAGS=-xHost
export CXXFLAGS=-xHost
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

rm data/*
###*********** RUNNING ************###
numactl --interleave=all ./bin/modeling-snap_end  --verbose --n1 512 --n2 512 --n3 512 --iter 504 --dshot 1 --first 1301 --last 1301
mv snapshot_final data/SB-final_snap_504.raw
numactl --interleave=all ./bin/modeling  --verbose --n1 512 --n2 512 --n3 512  --iter 504 --tb_thread_group_size 14 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 7 --tb_t_dim 3 --tb_num_wf 21 --dshot 1 --first 1301 --last 1301 -c 
mv TB-snapshot_final data/TB-final_snap_504.raw
rm lastshot_u1

numactl --interleave=all ./bin/modeling-snap_end  --verbose --n1 512 --n2 512 --n3 512 --iter 1000 --dshot 1 --first 1301 --last 1301
mv snapshot_final data/SB-final_snap_1000.raw
numactl --interleave=all ./bin/modeling  --verbose --n1 512 --n2 512 --n3 512  --iter 1000 --tb_thread_group_size 14 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 7 --tb_t_dim 3 --tb_num_wf 21 --dshot 1 --first 1301 --last 1301 -c
mv TB-snapshot_final data/TB-final_snap_1000.raw
rm lastshot_u1

numactl --interleave=all ./bin/modeling-snap_end  --verbose --n1 512 --n2 512 --n3 512 --iter 2000 --dshot 1 --first 1301 --last 1301
mv snapshot_final data/SB-final_snap_2000.raw
numactl --interleave=all ./bin/modeling  --verbose --n1 512 --n2 512 --n3 512  --iter 2000 --tb_thread_group_size 14 --tb_nb_thread_groups 4 --tb_th_x 1 --tb_th_y 2 --tb_th_z 7 --tb_t_dim 3 --tb_num_wf 21 --dshot 1 --first 1301 --last 1301 -c
mv TB-snapshot_final data/TB-final_snap_2000.raw
rm lastshot_u1

###*********** CHECKING RESULTS ***********###
./scripts/isnan data/TB-final_snap_504.raw
./scripts/isnan data/SB-final_snap_504.raw
./scripts/diff_to data/TB-final_snap_504.raw data/SB-final_snap_504.raw

./scripts/isnan data/TB-final_snap_1000.raw
./scripts/isnan data/SB-final_snap_1000.raw
./scripts/diff_to data/TB-final_snap_1000.raw data/SB-final_snap_1000.raw

./scripts/isnan data/TB-final_snap_2000.raw
./scripts/isnan data/SB-final_snap_2000.raw
./scripts/diff_to data/TB-final_snap_2000.raw data/SB-final_snap_2000.raw
