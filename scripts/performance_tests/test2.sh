#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=64:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=test2
#SBATCH --output=logs/test2_.%J.out
#SBATCH --error=logs/test2_.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** COMMENT *********###
# In this test we search for best SB parameters in terms of performance.
# 1)please change SLURM CONFIGURATION accordingly
#  partition,
# number of cpus-per-task should be equal to number of available cores.
# 2)please change OPENMP PARAMETERS accordingly
#OMP_NUM_THREADS should be equal to number of cpus
# 3)please change MODULES accordingly
#we need to load cmake,icpc modules from somewhere

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee test2.log) 2>&1
echo $hostname
#lscpu

###********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=48
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
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load icc/2020.2.254
module load cmake

##### Shot information #####
export shot=32896;# position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dx=10;

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_SB_1st=200
export NT_SB_2nd=200

##### Logs directory #####
mkdir ./logs
#########rm -rf ./logs/test2 #delete if existing
export logs_file="./logs/test2_attempt2.log"

##### File to modify #####
export file="./include/stencil/wave.h"

##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
    echo $i
    nx=${nx_arr[$i]}
    ny=${ny_arr[$i]}
    nz=${nz_arr[$i]}
    echo "grid nx=${nx}, ny=${ny}, nz=${nz}, th_z=${OMP_NUM_THREADS}"
#    for z in 1 `seq 2 10 64 `;do
    for x in 1 `seq 2 2 64 `;do     # 62 2 64
        for y in 1 `seq 2 2 64 `;do
            grid_str="${nx}_${ny}_${nz}_${x}_${y}"
            echo $grid_str;
            ##### Change cache blocking values #####
            sed -i "s/#define BLOCKX [0-9]\+/#define BLOCKX $x/" "$file"
            sed -i "s/#define BLOCKY [0-9]\+/#define BLOCKY $y/" "$file"
            sed -i "s/#define BLOCKZ [0-9]\+/#define BLOCKZ $z/" "$file"
            ##### COMPILATION #####
            mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
            CC=icc CXX=icpc cmake .
            make clean
            make VERBOSE=1
            make install
            #####################
            echo "Running SB"
            echo "Running 1st order"
            srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
            ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
            --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth \
             --order 1 --fmax $fmax \
            --dx $dx >> $logs_file

#            echo "Running 2nd order"
#            srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#            ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_2nd \
#            --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 2 --fmax $fmax \
#            --dx $dx >> $logs_path/log-SB_2nd-abc_$grid_str.log
        done
    done
#    done
done


