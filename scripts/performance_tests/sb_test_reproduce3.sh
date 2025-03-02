#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=reproduce_sb
#SBATCH --output=./logs/reproduce_sb.%J.out
#SBATCH --error=./logs/reproduce_sb.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # Don't use hyperthreading
########################################
###******** HORODATED LOG WRITING *********###
echo $HOSTNAME
lscpu

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
export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

###********** MODULES *********###
module load icc/2020.2.254
module load cmake
source /home/hltaief/pavel/intel/oneapi/setvars.sh

####################################################################
##### Shot information #####
export shot=32896
export src_depth=256
export fmax=8
export dx=10

###*********** Experiment setup ************###
nx=1024; ny=1024; nz=512
export NT_TB_2nd=530
export NT_SB_2nd=529
export NT_TB_1st=537
export NT_SB_1st=536
pwd

##### Logs directory #####
export logs_path="./logs/reproduce_sb"
mkdir -p "$logs_path"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Cache blocking to try #####
x=10;
y=22;
z=9999;

##### Run tests #####
grid_str="${nx}_${ny}_${nz}_${x}_${y}"
echo $grid_str

##### Change cache blocking values #####
sed -i "s/#define BLOCKX [0-9]\+/#define BLOCKX $x/" "$file"
sed -i "s/#define BLOCKY [0-9]\+/#define BLOCKY $y/" "$file"
sed -i "s/#define BLOCKZ [0-9]\+/#define BLOCKZ $z/" "$file"

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt
CC=icc CXX=icpc cmake .
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

#####################
#echo "Profiling with Callgrind..."
#callgrind_result_file="$logs_path/callgrind_$grid_str_$SLURM_JOB_ID.out"
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#valgrind --tool=callgrind --callgrind-out-file=$callgrind_result_file \
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
#--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax --dx $dx
#callgrind_annotate --auto=yes $callgrind_result_file > $logs_path/callgrind_report_$grid_str_$SLURM_JOB_ID.txt

echo "Job completed"
date