#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
####SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=reproduce_sb
#SBATCH --output=./logs/reproduce_sb.%J.out
#SBATCH --error=./logs/reproduce_sb.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading
########################################
###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee test2.log) 2>&1
echo $hostname
#lscpu

###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

# Set compiler flags
export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"  # -pg
export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

###********** MODULES *********###
# Load required modules
module load icc/2020.2.254
#source ~/.bashrc
#source /home/hltaief/pavel/intel/oneapi/setvars.sh
module load cmake

####################################################################
## Make sure the Intel environment is sourced properly
#source /home/hltaief/pavel/intel/oneapi/setvars.sh  # Sourcing Intel environment
## Set the correct compiler paths (ensure the Intel compiler is used)
#export CC=/home/hltaief/pavel/intel/oneapi/compiler/2024.1/bin/icx
#export CXX=/home/hltaief/pavel/intel/oneapi/compiler/2024.1/bin/icpx
####################################################################

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
CC=icc CXX=icpc cmake .
#cmake .
make clean
make VERBOSE=1
make install

#####################
echo "Running SB"
echo "Running 1st order"
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx >> $logs_path/log-SB_1st-abc_$grid_str.log

###### After execution, run gprof to generate the profiling report
gprof ./bin/modeling gmon.out > $logs_path/gprof_report.txt
