#!/bin/bash
###******** COMMENT *********###
# In this test we check performance with best SB,TB parameters for AMD MilanX.
# 1)please change SLURM CONFIGURATION accordingly
#  partition,
# number of cpus-per-task should be equal to number of available cores.
# 2)please change OPENMP PARAMETERS accordingly
#OMP_NUM_THREADS should be equal to number of cpus
# 3)please change MODULES accordingly
#we need to load cmake,icpc modules from somewhere

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee test.log) 2>&1
echo $hostname
#lscpu

###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=128
#export OMP_NUM_THREADS=64
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
#module load icc/2020.2.254
module load cmake

###********** Default SB, TB parameters *********###
export _CB_SIZE_X=16;
export _CB_SIZE_Y=4;
th_x_arr=(4 4 4)
th_y_arr=(2 2 2)
th_z_arr=(2 1 1)
num_wf_arr=(64 4 4)
tdim_arr=(7 7  7)

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=505
export NT_TB_2nd=502
export NT_SB_2nd=505

##### COMPILATION #####
###********** Set compiler flags *********###
export CFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export CXXFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export FFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
rm -rf ./logs/test1 #delete if existing
mkdir ./logs/test1
export logs_path=./logs/test1

##### Shot information #####
export shot=32896;# position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dx=10;

##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
#for i in $(seq 1 2); do
  echo $i
  nx=${nx_arr[$i]}
  ny=${ny_arr[$i]}
  nz=${nz_arr[$i]}
  echo "grid nx=${nx}, ny=${ny}, nz=${nz}, th_z=${OMP_NUM_THREADS}"
  grid_str="${nx}_${ny}_${nz}"
  th_x=${th_x_arr[$i]}
  th_y=${th_y_arr[$i]}
  th_z=${th_z_arr[$i]}
  t_dim=${tdim_arr[$i]}
  num_wf=${num_wf_arr[$i]}
  tgs=$((th_x * th_y*th_z))
  ###*********** SB ************###
  cbx=16; cby=4; cbz=9999;
#  echo "Running SB"
#  echo "Running 1st order"
#  ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
#  --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
#  --dx $dx --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/log-SB_1st-abc_$grid_str.log
#
#  echo "Running 2nd order"
#  ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_2nd \
#  --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 2 --fmax $fmax \
#  --dx $dx --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/log-SB_2nd-abc_$grid_str.log
  ###*********** TB ************##
  echo "Running TB"
  echo "Running 1st order"
  echo "num_th=${OMP_NUM_THREADS}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, t_dim=${t_dim}"
  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
  --src_depth $src_depth --order 1 --fmax $fmax --dx $dx >> $logs_path/log-TB_1st-abc_$grid_str.log;

  echo "Running 2nd order"
  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_2nd --tb_thread_group_size $tgs \
  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
  --src_depth $src_depth --order 2 --fmax $fmax --dx $dx >> $logs_path/log-TB_2nd-abc_$grid_str.log;
done
