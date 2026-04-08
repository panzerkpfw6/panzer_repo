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
#rm *snapshot*
#exit 0
rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
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
echo "Test 2. compare_wavefields_and_sismos_for_SB_TB"
export OMP_NUM_THREADS=4 #4
export TIME_TB_1st=505 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_SB_1st=505 #@pavel in SB the nt should one time less than in correponding TB.
### grid size 256*256*256
nx=256;ny=256;nz=256;
nx=200;ny=256;nz=160;

x=2; y=2; z=1; t=7; w=20; tgs=4;
export first=41100;
export first=41099;
export last=41100;
export src_depth=20;
export rcv_depth=4;

#mkdir ./data/sismos_sb
#rm ./data/sismos_sb/*
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2 --dshot 1 \
#  --first $shot --last $shot --src_depth $src_depth --rcv_depth $rcv_depth --drcv 1 --order 1 --fmax 8 --rec_sismos 1;
#mv ./data/sismos_${shot}.raw ./data/sismos_sb/sismos_${shot}.raw

  
mkdir ./data/sismos_tb
rm ./data/sismos_tb/*
./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
 --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $first --last $last -c \
 --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax 8 --rec_sismos 1;
mv ./data/sismos_${shot}.raw ./data/sismos_tb/sismos_${shot}.raw

./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505
./scripts_useful/diff_to ./data/sismos_sb/sismos_${shot}.raw ./data/sismos_tb/sismos_${shot}.raw
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