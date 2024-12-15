#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=test_default_pars
#SBATCH --output=logs/test1_.%J.out
#SBATCH --error=logs/test1_.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** COMMENT *********###
# In this test we check performance with best SB,TB parameters for AMD MilanX.
# 1)please change SLURM CONFIGURATION accordingly
#  partition,
# number of cpus-per-task should be equal to number of available cores.
# 2)please change OPENMP PARAMETERS accordingly
#OMP_NUM_THREADS should be equal to number of cpus
# 3)please change MODULES accordingly
#we need to load cmake,icpc modules from somewhere

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee test.log) 2>&1
echo $hostname
#lscpu

###********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=128
export OMP_NUM_THREADS=64
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
module load icc/2020.2.254
module load cmake

###********** Default SB, TB parameters *********###
export _CB_SIZE_X=8;export _CB_SIZE_Y=1;
th_x_arr=(8 4 4)
th_y_arr=(2 2 2)
th_z_arr=(1 1 1)
tdim_arr=(7 7  7)
num_wf_arr=(64 20 20)

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=504
export NT_TB_2nd=502
export NT_SB_2nd=501

##### COMPILATION #####
cd ./include/stencil
sed -i "s/#define BLOCKX [[:digit:]]\+[ ]*;/#define BLOCKX 8;/g" wave.h;


sed -i "s/#define BLOCKX [[:digit:]]\+[ ]*;/#define BLOCKX 8;/g" ./include/stencil/wave.h;

./include/stencil/wave.h
sed -i "s/#define BLOCKX [[:digit:]]\+[ ]*;/define BLOCKX $_CB_SIZE_X;/g" ./include/stencil/wave.h;

pwd
cd ./SB_abc;
sed -i "s/#define BLOCKX [[:digit:]]\+[ ]*;/define BLOCKX $_CB_SIZE_X;/g" SB_kernel.h;
sed -i "s/#define BLOCKY = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_Y = $_CB_SIZE_Y;/g" SB_kernel.h;
icpc -xHost -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_1st-abc.out

cd ../SB_order2_abc;
sed -i "s/const Myint _CB_SIZE_X = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_X = $_CB_SIZE_X  ;/g" SB_kernel.h;
sed -i "s/const Myint _CB_SIZE_Y = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_Y = $_CB_SIZE_Y;/g" SB_kernel.h;
icpc -xHost -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_2nd-abc.out

cd ../TB_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_TB_kernel.cpp -o ../TB_1st-abc.out
cd ../TB_order2_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2 -O3 -I. test_TB_kernel.cpp -o ../TB_2nd-abc.out
cd ..
#####################
##### Logs directory #####
mkdir ./logs
rm -rf ./logs/test1 #delete if existing
mkdir ./logs/test1
export logs_path=./logs/test1

##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
  echo $i
  nx=${nx_arr[$i]}
  ny=${ny_arr[$i]}
  nz=${nz_arr[$i]}
  echo "grid nx=${nx}, ny=${ny}, nz=${nz}, th_z=${OMP_NUM_THREADS}"
  grid_str="${nx}_${ny}_${nz}"
  th_x=${th_x_arr[$i]}
  th_y=${th_y_arr[$i]}
  th_z=${th_z_arr[$i]}
  t_dim=${tdim_arr[$i]}
  num_wf=${num_wf_arr[$i]}
  ###*********** SB ************###
  echo "Running SB"
  echo "Running 1st order"
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 ./SB_1st-abc.out  $nx $ny $nz $NT_SB_1st 3000 0 0 0 0 >> $logs_path/log-SB_1st-abc_$grid_str.log
  echo "Running 2nd order"
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 ./SB_2nd-abc.out  $nx $ny $nz $NT_SB_2nd 3000 0 0 >> $logs_path/log-SB_2nd-abc_$grid_str.log
  ###*********** TB ************##
  echo "Running TB"
  echo "Running 1st order"
  echo "num_th=${OMP_NUM_THREADS}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, t_dim=${t_dim}"
  srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --unbuffered numactl --interleave=all ./TB_1st-abc.out $nx $ny $nz $NT_TB_1st 3000 $th_x $th_y $th_z $t_dim $num_wf 0 0 0 0 >> "${logs_path}/log-TB_1st-abc_${grid_str}.log"
  echo "Running 2nd order"
  srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --unbuffered numactl --interleave=all ./TB_2nd-abc.out  $nx $ny $nz $NT_TB_2nd 3000 $th_x $th_y $th_z $t_dim $num_wf 0 0 >> "${logs_path}/log-TB_2nd-abc_${grid_str}.log"
done


