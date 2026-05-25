#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

# Set compiler flags
# export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
# export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
# export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

#export CFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export CXXFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export FFLAGS="-xHost -march=core-avx2 -mtune=core-avx2 -qopenmp -O3"


export CFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"
export CXXFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"
export FFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"

lscpu

###********** MODULES *********###
# Load required modules
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load icc/2020.2.254
# module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#source ~/.bashrc

module purge
module load intel/2025.3
module load cmake
# Source Intel oneAPI environment (this is REQUIRED for icx/icpx)
source /sw/rl9c/intel/oneapi/2025.3/setvars.sh --force
echo "Using compiler:"
which icx
which icpx
icx --version

# exit 0

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
export dh=10;
export dt=0.001;

###*********** Experiment setup ************###
nx=512; ny=512; nz=512;
export NT_TB_2nd=530  # @pavel in TB source injection starts from second time sample (Nothing happens for one dt). This is code feature.
export NT_SB_2nd=529  # @pavel in SB the nt should one time less than in corresponding TB.
export NT_TB_1st=537  # @pavel in TB source injection starts from second time sample (Nothing happens for one dt). This is code feature.
export NT_SB_1st=536  # @pavel in SB the nt should one time less than in corresponding TB.
export NT_SB_1st=100
pwd

##### Logs directory #####
mkdir logs
export logs_path="./logs/reproduce_sb"
mkdir "$logs_path"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Cache blocking to try #####
cbx=8; cby=6; cbz=9999;

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt  # Save the last CMakeCache.txt
#cmake .
CC=icx CXX=icpx cmake .
make clean
make VERBOSE=1
make install

##### Run tests #####
grid_str="${nx}_${ny}_${nz}_${cbx}_${cby}"
echo $grid_str;

#####################
echo "Running SB"
echo "Running 1st order"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dh --dy $dh --dz $dh --dt $dt --rec_sismos 0 \
--cbx $cbx --cby $cby --cbz $cbz;


# ###*********** SB ************###
# echo "Running SB"
# echo "Running 1st order"
# srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
# ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
# --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
# --dx $dh --dy $dh --dz $dh --dt $dt --rec_sismos 0 \
# --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/$logs_filename;

