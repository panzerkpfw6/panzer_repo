#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --cpus-per-task=128
#SBATCH --time=3:00:00
#SBATCH --partition=amd  #
#SBATCH --mem=200G
#SBATCH --job-name=test1_modeling
#SBATCH --output=logs/test1_modeling.%J.out
#SBATCH --error=logs/test1_modeling.%J.err
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** COMMENT *********###
# In this test we check performance with best SB,TB parameters for AMD MilanX.
# be aware that cache blocking should be different for 1st order and second order
# 1)please change SLURM CONFIGURATION accordingly
#  partition,
# number of cpus-per-task should be equal to number of available cores.
# 2)please change OPENMP PARAMETERS accordingly
#OMP_NUM_THREADS should be equal to number of cpus
# 3)please change MODULES accordingly
#we need to load cmake,icpc modules from somewhere

###******** srun command for job step execution of the script *********###
# srun --hint=nomultithread --nodes=1 --ntasks=1 --threads-per-core=1 --cpus-per-task=128 ./hello_world
###******** COMMENT *********###

###******** HORODATED LOG WRITING *********###
echo $hostname
lscpu

####********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=192
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
#export OMP_NESTED='True'
#export granularity=fine
#export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t
###********** Set compiler flags *********###
#export CFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export CXXFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"
#export FFLAGS="-march=core-avx2 -mtune=core-avx2 -qopenmp -O3"

######################################################
export OMP_NUM_THREADS=128;
export OMP_PLACES=cores;
export OMP_PROC_BIND=close;
export OMP_STACKSIZE=64M;
#export CFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
#export CXXFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
#export FFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"

# export CFLAGS="-march=znver4 -Ofast -fopenmp -fvectorize -floop-nest-optimize -floop-interchange -fno-math-errno -flto -mamdlibm"
# export CXXFLAGS="-march=znver4 -Ofast -fopenmp -fvectorize -floop-nest-optimize -floop-interchange -fno-math-errno -flto -mamdlibm"
# export FFLAGS="-march=znver4 -Ofast -fopenmp -fvectorize -floop-nest-optimize -floop-interchange -fno-math-errno -flto -mamdlibm"

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

##### Shot information #####
export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;
export dt=0.001;

###********** Default SB, TB parameters *********###
cbx_arr=(16  16 16)
cby_arr=(4  4  4)
cbz_arr=(9999  9999  9999)

### My recent parameter search
#th_x_arr_1st=(16 16 16)
#th_y_arr_1st=(2 2 2)
#th_z_arr_1st=(1 1 1)
#tdim_arr_1st=(7 7 7)
#num_wf_arr_1st=(192 192 192)

# PASC paper results
th_x_arr_1st=(8 4 4)
th_y_arr_1st=(2 2 2)
th_z_arr_1st=(1 1 1)
tdim_arr_1st=(7 7 7)
num_wf_arr_1st=(64 20 20)

### PASC-based results, my guess
#th_x_arr_1st=(8 2 2)
#th_y_arr_1st=(1 2 2)
#th_z_arr_1st=(1 1 1)
#tdim_arr_1st=(3 7 7)
#num_wf_arr_1st=(24 20 20)

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=505
#export NT_TB_1st=4001
#export NT_SB_1st=4001

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icx CXX=icpx cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
export logs_path=./logs/test1_modeling
export logs_filename="test1_fm_sbatch.log"
######rm -rf $logs_path
mkdir $logs_path

##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
# for i in $(seq 2 $len); do
# for i in $(seq 1 2); do
  echo $i
  nx=${nx_arr[$i]}
  ny=${ny_arr[$i]}
  nz=${nz_arr[$i]}
  cbx=${cbx_arr[$i]}
  cby=${cby_arr[$i]}
  cbz=${cbz_arr[$i]}
  echo "grid nx=${nx}, ny=${ny}, nz=${nz}, OMP_NUM_THREADS=${OMP_NUM_THREADS}"
  grid_str="${nx}_${ny}_${nz}"
  echo "cbx=${cbx}, cby=${cby}, cbz=${cbz}"

  ###*********** SB ************###
 echo "Running SB"
 echo "Running 1st order"
 srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
 ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
 --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
 --dx $dh --dy $dh --dz $dh --dt $dt --rec_sismos 0 \
 --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/$logs_filename;

  ###*********** TB ************##
  echo "Running TB"
  echo "Running 1st order"
  th_x=${th_x_arr_1st[$i]}
  th_y=${th_y_arr_1st[$i]}
  th_z=${th_z_arr_1st[$i]}
  t_dim=${tdim_arr_1st[$i]}
  num_wf=${num_wf_arr_1st[$i]}
  tgs=$((th_x * th_y*th_z))
  echo "num_th=${OMP_NUM_THREADS},th_x=${th_x},th_y=${th_y},th_z=${th_z},num_wf=${num_wf},t_dim=${t_dim},tgs=${tgs}"

#  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 numactl --interleave=all \
  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
  --src_depth $src_depth --order 1 --fmax $fmax --dx $dx --rec_sismos 0 >> $logs_path/$logs_filename;
done
