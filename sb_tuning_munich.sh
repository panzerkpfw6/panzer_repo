#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --nodelist=gbt04
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=tb_search
#SBATCH --output=./logs/tb_search.%J.out
#SBATCH --error=./logs/tb_search.%J.err
#SBATCH --hint=nomultithread
########################################
export TIME_TB_1st=505
export TIME_SB_1st=504
nt=500
mode=2
cluster_name="Milan"
th_x_arr=(1 2 4 8 16 32)
th_y_arr=(1 2)
th_z_arr=(1 2)
num_wf_arr=(2 4 8 12 16 20 24 32 64)
tdim_arr=(3 5 7 15) # suits for 512 domain size
########### Grid sizes search space
nx_arr=(  512  1024  2048  2048 )
ny_arr=(  512  1024  2048  2048 )
nz_arr=(  512  512   512   1024 )
#x_arr=(   512 )
#ny_arr=(   512 )
#nz_arr=(   512 )
###########
module load icc/2020.2.254

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
### COMPILE SB
echo "Compiling SB binaries"
cd ../SB; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_1st.out
cd ../SB_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_1st-abc.out
cd ../SB_order2; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_2nd.out
cd ../SB_order2_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_SB_kernel.cpp -o ../SB_2nd-abc.out
### COMPILE TB
echo "Compiling TB binaries"
echo "TB"
cd ../TB; icpc -qopenmp -march=core-avx2 -mtune=core-avx2 -O3 -I. test_TB_kernel.cpp -o ../TB_1st.out
echo "TB_abc"
cd ../TB_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2  -O3 -I. test_TB_kernel.cpp -o ../TB_1st-abc.out
echo "TB_order2"
cd ../TB_order2; icpc -qopenmp -march=core-avx2 -mtune=core-avx2 -O3 -I. test_TB_kernel.cpp -o ../TB_2nd.out
echo "TB_order2_abc"
cd ../TB_order2_abc; icpc -qopenmp -march=core-avx2 -mtune=core-avx2 -O3 -I. test_TB_kernel.cpp -o ../TB_2nd-abc.out
cd ..

###******** HORODATED LOG WRITING *********###
mkdir ./logs
mkdir ./logs/tb_param_search
export logs_path=./logs/tb_param_search/$SLURM_JOB_PARTITION
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
    for th_x in "${th_x_arr[@]}"; do
      for th_y in "${th_y_arr[@]}"; do
        for th_z in "${th_z_arr[@]}"; do
          tgs=$((th_x*th_y*th_z))
          if (( num_th < tgs )); then
            echo "SKIPPED num_th=${num_th}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, tdim=${tdim}"
            continue
          else
            for num_wf in "${num_wf_arr[@]}"; do
              for tdim in "${tdim_arr[@]}"; do
                echo "num_th=${num_th}, th_x=${th_x}, th_y=${th_y}, th_z=${th_z}, num_wf=${num_wf}, tdim=${tdim}"
                echo "${logs_path}/log-TB_1st_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}.log"
#                 srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./TB_1st-abc.out $nx $ny $nz $nt 300 $th_x $th_y $th_z $tdim $num_wf 0 0 0 0 >> "${logs_path}/log-TB_1st_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}_${grid_str}.log"
#                 srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./TB_2nd-abc.out $nx $ny $nz $nt 300 $th_x $th_y $th_z $tdim $num_wf 0 0 0 0 >> "${logs_path}/log-TB_2nd_${num_th}_${th_x}_${th_y}_${th_z}_${num_wf}_${tdim}_${grid_str}.log"
                srun --ntasks=1 --cpus-per-task=$num_th --hint=nomultithread numactl --interleave=all ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $nt --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $tdim --tb_num_wf $num_wf --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
              done
            done
          fi
        done
      done
    done
  done
done
echo "Halas"

