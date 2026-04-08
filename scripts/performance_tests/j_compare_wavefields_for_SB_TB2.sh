#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

# Set compiler flags
export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

###********** MODULES *********###
# Load required modules
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load icc/2020.2.254
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#source ~/.bashrc
module load cmake

###***********************###
# Make sure the Intel environment is sourced properly
# Set the correct compiler paths (ensure the Intel compiler is used)
#export CC=/opt/intel/oneapi/compiler/2023.2.4/linux/bin/icx
#export CXX=/opt/intel/oneapi/compiler/2023.2.4/linux/bin/icpx
###***********************###

##### Shot information #####
export shot=32896;  # position of the source in x,y coordinates. Check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dx=10;

###*********** Experiment setup ************###
nx=512; ny=512; nz=512;
export NT_TB_2nd=530  # @pavel in TB source injection starts from second time sample (Nothing happens for one dt). This is code feature.
export NT_SB_2nd=529  # @pavel in SB the nt should one time less than in corresponding TB.
export NT_TB_1st=537  # @pavel in TB source injection starts from second time sample (Nothing happens for one dt). This is code feature.
export NT_SB_1st=536  # @pavel in SB the nt should one time less than in corresponding TB.
export NT_SB_1st=100
pwd

##### Logs directory #####
export logs_path="./logs/reproduce_sb"
mkdir "$logs_path"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Cache blocking to try #####
x=10; y=22;

##### Run tests #####
grid_str="${nx}_${ny}_${nz}_${x}_${y}"
echo $grid_str;

##### Change cache blocking values #####
sed -i "s/#define BLOCKX [0-9]\+/#define BLOCKX $x/" "$file"
sed -i "s/#define BLOCKY [0-9]\+/#define BLOCKY $y/" "$file"

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt  # Save the last CMakeCache.txt
#cmake .
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

#####################
echo "Running SB"
echo "Running 1st order"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx
# >> $logs_path/log-SB_1st-abc_$grid_str.log
## After execution, run gprof to generate the profiling report
#gprof ./bin/modeling gmon.out > $logs_path/gprof_report.txt

#####################
#echo "Profiling with perf..."
#perf record -g ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
#--dx $dx
#perf report > $logs_path/perf_report.txt

#####################
#echo "Running memory check with valgrind..."
#valgrind --tool=memcheck --leak-check=full --show-reachable=yes --track-origins=yes \
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
#--dx $dx 2> $logs_path/valgrind_memcheck.txt