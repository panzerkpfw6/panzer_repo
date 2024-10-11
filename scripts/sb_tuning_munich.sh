#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --nodelist=gbt07
#SBATCH --time=72:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=sb_search
#SBATCH --output=./logs/sb_search.%J.out
#SBATCH --error=./logs/sb_search.%J.err
#SBATCH --hint=nomultithread
########################################
nt=504
mode=2
cluster_name="Milan"
########### Grid sizes search space
nx_arr=(  512  1024  2048  2048 )
ny_arr=(  512  1024  2048  2048 )
nz_arr=(  512  512   512   1024 )
x_arr=(   512 )
ny_arr=(   512 )
nz_arr=(   512 )
###########
module load icc/2020.2.254
module load cmake

if [ "$SLURM_JOB_PARTITION" = "7763" ]; then
    num_th_arr=(1 32 64 96 128)
#    num_th_arr=(128)
elif [ "$SLURM_JOB_PARTITION" = "7773X" ]; then
#    num_th_arr=(1 32 64 96 128)
    num_th_arr=(128)
else
    exit 1
fi

echo $HOSTNAME
echo $SLURM_JOB_PARTITION
module list
lscpu

###*********** compilation ************###
#OpenMP settings:
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t
module load icc/2020.2.254
module list
###******** HORODATED LOG WRITING *********###
mkdir ./logs
mkdir ./logs/sb_param_search
export logs_path=./logs/sb_param_search/$SLURM_JOB_PARTITION
mkdir $logs_path
exec > >(while read line; do echo "$(date): $line"; done | tee ./logs/tb_search_$SLURM_JOB_PARTITION.log) 2>&1
###*********** RUNNING SIMWAVE ************###
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
    ## in the code only BLOCKY,BLOCKZ are used. set up cache blocking.
    for x in 1 `seq 2 2 64 `;do
        for y in 1 `seq 2 2 64 `;do
          echo $x,$y;
#          sed -i "s/define BLOCKY [[:digit:]]\+[ ]*;/define BLOCKY $x;/g" ./include/stencil/wave.h;
#          sed -i "s/define BLOCKZ [[:digit:]]\+[ ]*;/define BLOCKZ $y;/g" ./include/stencil/wave.h;
          sed -i "s/.*#define BLOCKZ [[:digit:]].*/#define BLOCKZ $x/g" ./include/simwave/wave.h;
          sed -i "s/.*#define BLOCKY [[:digit:]].*/#define BLOCKY $x/g" ./include/simwave/wave.h;

          echo "compilation"
          mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
          CC=icc CXX=icpc cmake .
          make clean
          make VERBOSE=1
          make install

          echo "${logs_path}/log-SB_1st_${x}_${y}.log"
          srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --mode $mode  --dshot 1 --first 1301 --last 1301 >> "${logs_path}/log-SB_1st_${x}_${y}.log"
        done
    done
  done
done
echo "Halas"

