#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=48:00:00
#SBATCH --nodelist=gbt07
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=test3
#SBATCH --output=logs/test3_.%J.out
#SBATCH --error=logs/test3_.%J.err
#SBATCH --cpus-per-task=128
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
exec > >(while read line; do echo "$(date): $line"; done | tee test3.log) 2>&1
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

###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load icc/2020.2.254
module load cmake

###********** TB parameters diapason *********###
num_th_arr=(128)
#th_x_arr=(1 2 4 8 16 32)
th_x_arr=(8 16 32)
th_y_arr=(1 2 4 8 16 32)
th_z_arr=(1 2 4 8 16 32)
num_wf_arr=(2 4 8 12 16 20 24 32 64)
tdim_arr=(3 5 7 15) # suits for 512 domain size

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=504
export NT_TB_2nd=502
export NT_SB_2nd=501

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
#####rm -rf ./logs/test3 #delete if existing
mkdir ./logs/test3
export logs_path=./logs/test3

##### Shot information #####
export shot=32896;# position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dx=10;

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
                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread --threads-per-core=1 numactl --interleave=all \
                ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
                --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
                --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
                --src_depth $src_depth --order 1 --fmax $fmax --dx $dx >> "${logs_path}/log-TB_1st_abc_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${t_dim}_${grid_str}.log";

                echo "Running 2nd order"
                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread --threads-per-core=1 numactl --interleave=all \
                ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_2nd --tb_thread_group_size $tgs \
                --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
                --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
                --src_depth $src_depth --order 2 --fmax $fmax --dx $dx >> "${logs_path}/log-TB_2nd_abc_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${t_dim}_${grid_str}.log";

#                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread --threads-per-core=1 numactl --interleave=all ./TB_1st-abc.out $nx $ny $nz $NT_TB_1st 3000 $th_x $th_y $th_z $tdim $num_wf 0 0 0 0 >> "${logs_path}/log-TB_1st_abc_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}_${grid_str}.log"
#                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread --threads-per-core=1 numactl --interleave=all ./TB_2nd-abc.out $nx $ny $nz $NT_TB_2nd 3000 $th_x $th_y $th_z $tdim $num_wf 0 0 >> "${logs_path}/log-TB_2nd_abc_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}_${grid_str}.log"
              done
            done
          fi
        done
      done
    done
  done
done
