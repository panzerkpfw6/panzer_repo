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

###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load icc/2020.2.254
module load cmake

###********** SB, TB parameters diapason *********###

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=504
export NT_TB_2nd=502
export NT_SB_2nd=501

##### Logs directory #####
mkdir ./logs
rm -rf ./logs/test2 #delete if existing
mkdir ./logs/test2
export logs_path=./logs/test2

##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
    echo $i
    nx=${nx_arr[$i]}
    ny=${ny_arr[$i]}
    nz=${nz_arr[$i]}
    echo "grid nx=${nx}, ny=${ny}, nz=${nz}, th_z=${OMP_NUM_THREADS}"
    for x in 1 `seq 2 2 64 `;do     # 62 2 64
        for y in 1 `seq 2 2 64 `;do
            grid_str="${nx}_${ny}_${nz}_${x}_${y}"
            echo $grid_str;
            ##### COMPILATION #####
            cd ./SB_abc;
            sed -i "s/const Myint _CB_SIZE_X = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_X = $x  ;/g" SB_kernel.h;
            sed -i "s/const Myint _CB_SIZE_Y = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_Y = $y;/g" SB_kernel.h;
            icpc -xHost -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_1st-abc.out

            cd ../SB_order2_abc;
            sed -i "s/const Myint _CB_SIZE_X = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_X = $x  ;/g" SB_kernel.h;
            sed -i "s/const Myint _CB_SIZE_Y = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_Y = $y;/g" SB_kernel.h;
            icpc -xHost -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_2nd-abc.out
            cd ..
            #####################
            echo "Running 1st order"
            srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 ./SB_1st-abc.out  $nx $ny $nz $NT_SB_1st 3000 0 0 0 0 >> $logs_path/log-SB_1st-abc_$grid_str.log
            echo "Running 2nd order"
            srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 ./SB_2nd-abc.out  $nx $ny $nz $NT_SB_2nd 3000 0 0 >> $logs_path/log-SB_2nd-abc_$grid_str.log
        done
    done
done


