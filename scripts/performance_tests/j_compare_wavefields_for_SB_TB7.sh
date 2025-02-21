#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

###********** MODULES *********###
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
module load cmake
source /media/plotnips/sdd1/soft/intel/oneapi/setvars.sh

##### Shot information #####
export shot=32896
export src_depth=256
export fmax=8
export dx=10

###*********** Experiment setup ************###
nx=512; ny=512; nz=512
export NT_SB_1st=100
cd /media/plotnips/sdd1/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/stencil_rtm_ordering
pwd

##### Logs directory #####
export logs_path="./logs/reproduce_sb"
mkdir -p "$logs_path"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Cache blocking to try #####
x=10; y=22

##### Run tests #####
grid_str="${nx}_${ny}_${nz}_${x}_${y}"
echo $grid_str

##### Change cache blocking values #####
sed -i "s/#define BLOCKX [0-9]\+/#define BLOCKX $x/" "$file"
sed -i "s/#define BLOCKY [0-9]\+/#define BLOCKY $y/" "$file"

##### COMPILATION #####
export CFLAGS="-O3 -qopenmp -g"
export CXXFLAGS="-O3 -qopenmp -g"
export FFLAGS="-O3 -qopenmp -g"
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

#####################
echo "Running SB"
echo "Running 1st order"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx >> $logs_path/log-SB_1st-abc_3_$grid_str.log

#####################
echo "Profiling with VTune..."
vtune -collect hotspots -r $logs_path/vtune_hotspots_$grid_str ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax --dx $dx