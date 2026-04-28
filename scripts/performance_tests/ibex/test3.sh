#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=workq  #
#SBATCH --job-name=find_pars_TB
#SBATCH --output=logs/test3_.%J.out
#SBATCH --error=logs/test3_.%J.err
#SBATCH --cpus-per-task=192
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** COMMENT *********###
# In this test we search for best TB parameters in terms of performance.
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
export OMP_PLACES=cores;
export OMP_PROC_BIND=close;
export OMP_STACKSIZE=64M;
export OMP_NUM_THREADS=192;
export CFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export CXXFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export FFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"

###********** MODULES *********###
module load intel-oneapi/2023.1.0
module load cmake

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

###********** TB parameters diapason *********###
num_th_arr=(192)
th_x_arr=(1 2 4 8 16 32)
th_y_arr=(1 2 4 8 16 32)
th_z_arr=(1 2 4 8 16 32)
num_wf_arr=(2 4 8 12 16 20 24 32 64 128 192)
tdim_arr=(3 5 7 15) # suits for 512 domain size

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=200
export NT_TB_2nd=200

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs

######### create log
export logs_file="./logs/test3_modeling.log"
rm $logs_file
lscpu >> logs_file
echo $hostname

##### Shot information #####
export shot=32896;# position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;

##### Run tests #####
# iterate on parameters values
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
  nx=${nx_arr[$i]}
  ny=${ny_arr[$i]}
  nz=${nz_arr[$i]}
  echo "grid nx=${nx}, ny=${ny}, nz=${nz}, num_th=${OMP_NUM_THREADS}"
  grid_str="${nx}_${ny}_${nz}"
  for num_th in "${num_th_arr[@]}"; do
    export OMP_NUM_THREADS=$num_th
    for th_x in "${th_x_arr[@]}"; do
      for th_y in "${th_y_arr[@]}"; do
        for th_z in "${th_z_arr[@]}"; do
          tgs=$((th_x*th_y*th_z))
          if (( num_th < tgs )); then
            echo "SKIPPED num_th=${num_th}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, tdim=${t_dim}"
            continue
          else
            for num_wf in "${num_wf_arr[@]}"; do
              for t_dim in "${tdim_arr[@]}"; do
                echo "num_th=${num_th}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, t_dim=${t_dim}"
                echo "${logs_path}/log-TB_1st_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${t_dim}.log"

                echo "Running 1st order"
                echo "num_th=${OMP_NUM_THREADS}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, t_dim=${t_dim}"
                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
                ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
                --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
                --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
                --src_depth $src_depth --order 1 --fmax $fmax --dx $dh --dy $dh --dz $dh >> $logs_file

              done
            done
          fi
        done
      done
    done
  done
done
