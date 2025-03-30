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
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#module load icc/2020.2.254
module load cmake
##### COMPILATION #####
rm ./bin/modeling
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

############################### parameters
nx=256;ny=256;nz=256;
#nx=512;ny=512;nz=512;
export TIME_SB_1st=537; #@pavel in SB the nt should one time less than in correponding TB.
export TIME_SB_1st=100;
export TIME_TB_1st=100;
export TIME_TB_1st=537;
export TIME_TB_2nd=537;

export TIME_SB_1st=80;
export TIME_SB_2nd=80;


dt=0.001;
shot=16447;  # position of the source in x,y coordinates.check ./data/acquisition.txt
shot=32638;  # position of the source in x,y coordinates.check ./data/acquisition.txt
shot=32896;
#shot=100;  # position of the source in x,y coordinates.check ./data/acquisition.txt
src_depth=128;
x=2; y=2; z=1; t=7; w=20; tgs=4;
export OMP_NUM_THREADS=12;
############################### small test
export TIME_SB_1st=30;
#nx=1024;ny=1024;nz=512;
############################### compare_wavefields_for_SB_TB_2nd_order
## SB
#gdb --args ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_SB_2nd \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth \
#--cbx 2 --cby 2 --cbz 2 --order 2 --fmax 8 --dx 10;

## works
./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_SB_2nd \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth \
--cbx 8 --cby 1 --cbz 9999 --order 2 --fmax 8 --dx 10;
#
## works
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_2nd --tb_thread_group_size $tgs \
# --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
# --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 2 --fmax 8;
#
##gdb --args ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_2nd --tb_thread_group_size $tgs \
##  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
##  --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 2 --fmax 8;
#
# ./scripts_useful/diff_to ./snapshot_TB2nd_530 ./snapshot_SB2nd_530;
 ############################### compare_wavefields_for_SB_TB_1st_order
 # works
./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_SB_1st \
 --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth \
--cbx 8 --cby 1 --cbz 9999  --order 1 --fmax 8 --dx 10;

### debugging with gdb
#gdb --args ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_SB_1st \
# --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax 8 --dx 10;

# works
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
#  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#  --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 1 --fmax 8 --dx 10;

#gdb --args ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
#  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#  --tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c --src_depth $src_depth --order 1 --fmax 8 --dx 10;
##
#./scripts_useful/diff_to ./snapshot_TB1st_537 ./snapshot_SB1st_537;

############################################################
