#!/bin/bash

###******** HORODATED LOG WRITING *********###
echo $hostname
lscpu
nvidia-smi
########### CPU OPTIONS ####################
export OMP_NUM_THREADS=128
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_WAIT_POLICY=active
export KMP_BLOCKTIME=0
export KMP_AFFINITY=granularity=core,compact=1,1
export OMP_MAX_ACTIVE_LEVELS=1
export OMP_DYNAMIC=false

# export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3 -DUSE_CUDA"
# export CXXFLAGS="$CFLAGS"
# export FFLAGS="$CFLAGS"


########### GPU OPTIONS ####################
unset CFLAGS
unset CXXFLAGS
unset FFLAGS
unset CUDAFLAGS
unset NVCCFLAGS

######################################################
######################################################
###********** MODULES *********###
module purge
module load cuda/12.4.1
module load intel/2022.3
# source /sw/rl9c/intel/2025/compiler/2025.3/env/vars.sh --force   # or the correct path
module load cmake

echo "Using compiler:"
which icx
which icpx
icx --version

##### COMPILATION #####
# CPU version
# mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
## CC=icx CXX=icpx cmake .     # regular cpu

# GPU version
rm -rf CMakeFiles CMakeCache.txt
cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER="$(which icx)" \
  -DCMAKE_CXX_COMPILER="$(which icpx)" \
  -DCMAKE_CUDA_COMPILER="$(which nvcc)" \
  -DCMAKE_CUDA_ARCHITECTURES=70 \
  -DUSE_CUDA=ON \
  .
make clean
make VERBOSE=1
make install

# we need also change this in CMakeLists.txt:
# project(panzer_repo LANGUAGES C CXX CUDA)
exit 1

##### Shot information #####
export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;
export dt=0.001;

###********** Default SB, TB parameters *********###
cbx_arr=(16  16 16)
cby_arr=(4  4  4)
cbz_arr=(9999  9999  9999)

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=505
#export NT_TB_1st=4001
#export NT_SB_1st=4001

###############################
#echo "test_SB"
#nx=128;ny=256;nz=512;
#nt=10;  dt=0.001;
#export shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
#export src_depth=256;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --mode 2 --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --mode 2  --dshot 1 --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;

###############################
#echo "test_TB"
#nx=128;ny=256;nz=512;
#nt=57; dt=0.001;
#x=2; y=2; z=1; t=7; w=20; tgs=4;
#export OMP_NUM_THREADS=4
#export shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
#export src_depth=256;
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 1 --fmax 8;
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 2 --fmax 8;

################################
#echo "compare_wavefields_for_SB_TB"
#export OMP_NUM_THREADS=4 #4
#export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
#export TIME_SB_1st=505 #@pavel in SB the nt should one time less than in correponding TB.
#### grid size 256*256*256
#nx=256;ny=256;nz=256;
#x=2; y=2; z=1; t=7; w=20; tgs=4;
#export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
#export src_depth=128;
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
# --src_depth $src_depth --order 1 --fmax 8;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
#  --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
#./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505

################################
# echo "Test 2. compare_wavefields_and_sismos_for_SB_TB"
# export OMP_NUM_THREADS=4 #4
# export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
# export TIME_SB_1st=505 #@pavel in SB the nt should one time less than in correponding TB.
# ### grid size 256*256*256
# nx=256;ny=256;nz=256;
# nx=200;ny=256;nz=160;

# x=2; y=2; z=1; t=7; w=20; tgs=4;
# export first=41100;
# export first=41099;
# export last=41100;
# export src_depth=20;
# export rcv_depth=4;

# #mkdir ./data/sismos_sb
# #rm ./data/sismos_sb/*
# #./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
# #  --first $shot --last $shot --src_depth $src_depth --rcv_depth $rcv_depth --drcv 1 --order 1 --fmax 8 --rec_sismos 1;
# #mv ./data/sismos_${shot}.raw ./data/sismos_sb/sismos_${shot}.raw

# exit 1
  
# mkdir ./data/sismos_tb
# rm ./data/sismos_tb/*
# ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
#  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#  --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $first --last $last -c \
#  --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax 8 --rec_sismos 1;
# mv ./data/sismos_${shot}.raw ./data/sismos_tb/sismos_${shot}.raw

# ./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505
# ./scripts_useful/diff_to ./data/sismos_sb/sismos_${shot}.raw ./data/sismos_tb/sismos_${shot}.raw
#
################################
#echo "compare_wavefields_for_SB_TB_2nd_order"
#export OMP_NUM_THREADS=4
#export TIME_TB_2nd=514 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
#export TIME_SB_2nd=514 #@pavel in SB the nt should one time less than in correponding TB.
#### grid size 256*256*256
#nx=256;ny=256;nz=256;
#x=2; y=2; z=1; t=7; w=20; tgs=4;
#export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
#export src_depth=128;
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_2nd --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
# --src_depth $src_depth --order 2 --fmax 8;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_2nd --mode 2 --dshot 1 \
#  --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
#./scripts_useful/diff_to ./snapshot_TB2nd_514 ./snapshot_SB2nd_514

###############################
#echo "compare_wavefields_for_SB_1st_2nd_order"
#export OMP_NUM_THREADS=4
#export TIME_SB_1st=520 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
#export TIME_SB_2nd=520 #@pavel in SB the nt should one time less than in correponding TB.
#### grid size 256*256*256
#nx=256;ny=256;nz=256;
#export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
#export src_depth=128;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
#  --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 1 --fmax 8;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_2nd --mode 2 --dshot 1 \
#  --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
#./scripts_useful/diff_to ./snapshot_SB1st_520 ./snapshot_SB2nd_520

################################
#echo "test_sismos_options_for_SB"
################################
#echo "test_sismos_options_for_TB"
################################
#echo "compare_sismos_for_SBabc_TBabc"
################################
#echo "test_sismos_options_for_SB"
##############################