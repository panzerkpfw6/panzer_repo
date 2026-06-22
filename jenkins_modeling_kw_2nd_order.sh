#!/bin/bash
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
module load cmake

##### COMPILATION #####
pwd
rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved

##### build the project
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

#####
echo "compare_wavefields_for_SB_TB_2nd_order"
export OMP_NUM_THREADS=4
export TIME_TB_2nd=514 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_SB_2nd=514 #@pavel in SB the nt should one time less than in correponding TB.
### grid size 256*256*256
nx=256;ny=256;nz=256;
x=2; y=2; z=1; t=7; w=20; tgs=4;
export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=128;

echo "TB"
./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_2nd --tb_thread_group_size $tgs \
 --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
 --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
 --src_depth $src_depth --order 2 --fmax 8;

 echo "SB"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_2nd --mode 2 --dshot 1 \
  --first $shot --last $shot --src_depth $src_depth --drcv 1 --order 2 --fmax 8;
./scripts_useful/diff_to ./snapshot_TB2nd_514 ./snapshot_SB2nd_514