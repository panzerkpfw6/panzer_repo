#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=64  # Optimized for physical cores (assuming 64 on Milan-X)
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_MAX_ACTIVE_LEVELS=1  # Replaces deprecated OMP_NESTED
export granularity=fine
export KMP_AFFINITY=scatter  # Better core utilization
export KMP_HW_SUBSET=1t

###********** MODULES *********###
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
module load cmake
source /media/plotnips/sdd1/soft/intel/oneapi/setvars.sh

##### Shot information #####
#export shot=32896
export shot=524800
export src_depth=256
export fmax=8
export dx=10

###*********** Experiment setup ************###
nx=512; ny=512; nz=512
nx=1024; ny=1024; nz=512
grid_str="${nx}_${ny}_${nz}_${x}_${y}"
echo $grid_str
export NT_SB_1st=100
cd /media/plotnips/sdd1/Dropbox/PhD_proposal/work_with_david/Exawave_3_handover/stencil_rtm_ordering
pwd

##### Logs directory #####
export logs_path="./logs/reproduce_sb"
mkdir -p "$logs_path"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Cache blocking to try #####
x=10; y=22; z=9999;  # Optimized for cache and threads
echo "Updating cache blocking: BLOCKX=$x, BLOCKY=$y, BLOCKZ=$z"
sed -i "s/#define BLOCKX [0-9]\+/#define BLOCKX $x/" "$file"
sed -i "s/#define BLOCKY [0-9]\+/#define BLOCKY $y/" "$file"
sed -i "s/#define BLOCKZ [0-9]\+/#define BLOCKZ $z/" "$file"


##### COMPILATION  (performance) #####
export CFLAGS="-O3 -qopenmp -g"
export CXXFLAGS="-O3 -qopenmp -g"
export FFLAGS="-O3 -qopenmp -g"

##### COMPILATION  (better debugging) #####
#export CFLAGS="-O0 -qopenmp -g"  # Changed -O3 to -O0 for better debugging
#export CXXFLAGS="-O0 -qopenmp -g"
#export FFLAGS="-O0 -qopenmp -g"

CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

#####################
echo "Running SB"
echo "Running 1st order"
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
#--dx $dx >> $logs_path/log-SB_1st-abc_5_$grid_str.log

./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx

#gdb --args ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
#--dx $dx



#####################
#echo "Profiling with VTune..."
#vtune_result_dir="$logs_path/vtune_hotspots_$grid_str_$(date +%s)"  # Unique timestamp
#/media/plotnips/sdd1/soft/intel/oneapi/vtune/2024.1/bin64/vtune -collect hotspots -r "$vtune_result_dir" \
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax --dx $dx
#/media/plotnips/sdd1/soft/intel/oneapi/vtune/2024.1/bin64/vtune -report hotspots -r "$vtune_result_dir" > "$logs_path/vtune_report_$grid_str_$(date +%s).txt"

#