#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=128
#SBATCH --time=24:00:00
#SBATCH --partition=amd  #
#SBATCH --mem=200G
#SBATCH --job-name=find_pars_SB_rome
#SBATCH --output=logs/test2_.%J.out
#SBATCH --error=logs/test2_.%J.err
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
echo $hostname
lscpu

###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=128;
export OMP_PLACES=cores;
export OMP_PROC_BIND=close;
export OMP_STACKSIZE=64M;
export OMP_NUM_THREADS=192;
export CFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"
export CXXFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"
export FFLAGS="-march=znver2 -m64 -Ofast -ffast-math -qopenmp -O3"

###********** MODULES *********###
module purge
module load intel/2025.3
module load cmake
# Source Intel oneAPI environment (this is REQUIRED for icx/icpx)
source /sw/rl9c/intel/oneapi/2025.3/setvars.sh --force
echo "Using compiler:"
which icx
which icpx
icx --version

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icx CXX=icpx cmake .
make clean
make VERBOSE=1
make install

##### Shot information #####
export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;
export dt=0.001;

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_SB_1st=200
export NT_SB_2nd=200

##### Logs directory #####
mkdir ./logs
export logs_file="./logs/test2_modeling.log"
# rm $logs_file
lscpu >> logs_file
echo $hostname

##### Shot information #####
export shot=32896;# position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;

##### Run tests #####
len=${#nx_arr[@]}
export cbz=9999;
for i in $(seq 0 $len); do
    echo $i
    nx=${nx_arr[$i]}
    ny=${ny_arr[$i]}
    nz=${nz_arr[$i]}
    echo "grid nx=${nx}, ny=${ny}, nz=${nz}, th_z=${OMP_NUM_THREADS}"
#    for z in 1 `seq 2 10 64 `;do
    for cbx in 1 `seq 2 2 64 `;do     # 62 2 64
        for cby in 1 `seq 2 2 64 `;do
            grid_str="${nx}_${ny}_${nz}_${cbx}_${cby}"
            echo $grid_str;
            #####################
            echo "Running SB"
            echo "Running 1st order"
            srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
            ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
            --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth \
             --order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz \
            --dx $dh --dy $dh --dz $dh >> $logs_file
        done
    done
#    done
done


