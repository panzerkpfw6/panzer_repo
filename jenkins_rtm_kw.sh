#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=44
#export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t
###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#module load icc/2020.2.254
module load cmake
##### COMPILATION #####
pwd


rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
#rm ./data/*ilm*
#rm ./data/*img*
#rm ./data/*sismos*
#rm ./data/*snap*
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved

##### build the project with debugging symbols
#CC=icc CXX=icpc -DCMAKE_BUILD_TYPE=Debug cmake .
#CC=icc CXX=icpc -O0 -g3 -DCMAKE_BUILD_TYPE=Debug cmake .

##### build the project
CC=icc CXX=icpc cmake .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg" -DCMAKE_CXX_FLAGS="-pg" .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg -O2" -DCMAKE_CXX_FLAGS="-pg -O2" .

make clean
make VERBOSE=1
make install

################################
echo "Test 1. compare_rtm_result_for_SB_TB"
export OMP_NUM_THREADS=4 #4
export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_SB_1st=505 #@pavel in SB the nt should one time less than in correponding TB.
dt=0.001;fmax=11;
### grid size 256*256*256
nx=256;ny=256;nz=256;dh=25;
nx=200;ny=256;nz=160;dh=25;

x=2; y=2; z=1; t=7; w=20; tgs=4;
cbx=64;cby=22;cbz=9999;

first=1626; # position of the source:isx=160,isy=140
last=1627; # position of the source:isx=160,isy=140
export src_depth=5;
export rcv_depth=8;
##################### SB RTM workflow	#####################
#mkdir ./data/rtm_sb
#rm ./data/rtm_sb/*

./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz;

#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
#--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
#--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
#--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz;
#
#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
#--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
#--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
#--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz;

#mv ./data/sismos_${shot}.raw ./data/rtm_sb/sismos_${shot}.raw 
#########################################################

  
#mkdir ./data/rtm_tb
#rm ./data/rtm_tb/*
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
# --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax 8;
#mv ./data/sismos_${shot}.raw ./data/rtm_tb/sismos_${shot}.raw

#./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505
#./scripts_useful/diff_to ./data/rtm_sb/sismos_${shot}.raw ./data/rtm_tb/sismos_${shot}.raw
#
