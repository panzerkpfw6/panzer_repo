#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=workq
#SBATCH --job-name=test_default_pars
#SBATCH --output=logs/test1_.%J.out
#SBATCH --error=logs/test1_.%J.err
#SBATCH --cpus-per-task=192
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

###******** HORODATED LOG WRITING *********###
echo $hostname
lscpu

###********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=192
##export OMP_NUM_THREADS=64
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
export OMP_PLACES=cores;
export OMP_PROC_BIND=close;
export OMP_STACKSIZE=64M;
export OMP_NUM_THREADS=192;
export CFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export CXXFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"
export FFLAGS="-march=znver4 -dynamic -m64 -Ofast -ffast-math -fopenmp -O3"

###********** MODULES *********###
########module load intel/2024.2.1
module load intel-oneapi/2023.1.0
module load cmake

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
rm -rf ./logs/test1 #delete if existing
mkdir ./logs/test1
export logs_path=./logs/test1

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
timesteps=2200
nx=676;ny=676;nz=201;
#### Profile x=310. Salt3D. no mistake
first=20957;last=21023;
#### Profile x=?. Salt3D. no mistake
first=1;last=110;
dshot=450;
fmax=11;

./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
   >> $logs_path/log_model_water.log;
exit 1;

echo "Model data for RTM. water halfspace."
first=20468;last=20507;
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
   >> $logs_path/log_model_water.log;
 
#echo "Do Python filtering of real data." 
# 
#echo "Model data for RTM. salt3d."
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
#--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
# >> $logs_path/log_model_salt3d.log;
#
#echo "Perform RTM"
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
#--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
#>> $logs_path/log_rtm_salt3d.log;
#	
#
#echo "Gather images"
#srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 \
# --first $first --last $last -c --fwd_steps 3 --order 1 --src_depth 5 --rcv_depth 8 --drcv 1 --dir "./data" \
# >> $logs_path/log_gather_salt3d.log;





##### Run tests #####
len=${#nx_arr[@]}
for i in $(seq 0 $len); do
#for i in $(seq 0 1); do
  echo $i
  nx=${nx_arr[$i]}
  ny=${ny_arr[$i]}
  nz=${nz_arr[$i]}
  echo "grid nx=${nx}, ny=${ny}, nz=${nz}, OMP_NUM_THREADS=${OMP_NUM_THREADS}"
  grid_str="${nx}_${ny}_${nz}"
  ###*********** SB ************###
  echo "Running SB"
  echo "Running 1st order"
  cbx=64; cby=22; cbz=9999;
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
  ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
  --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
  --dx $dx --cbx $cbx --cby $cby --cbz $cbz  >> $logs_path/log-SB_1st-abc_$grid_str.log

  echo "Running 2nd order"
#  cbx=16; cby=4; cbz=9999;
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
  ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_2nd \
  --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 2 --fmax $fmax \
  --dx $dx --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/log-SB_2nd-abc_$grid_str.log
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

  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
  --src_depth $src_depth --order 1 --fmax $fmax --dx $dx >> $logs_path/log-TB_1st-abc_$grid_str.log;
#
  echo "Running 2nd order"
  th_x=${th_x_arr_2nd[$i]}
  th_y=${th_y_arr_2nd[$i]}
  th_z=${th_z_arr_2nd[$i]}
  t_dim=${tdim_arr_2nd[$i]}
  num_wf=${num_wf_arr_2nd[$i]}
  tgs=$((th_x*th_y*th_z))
  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_2nd --tb_thread_group_size $tgs \
  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
  --src_depth $src_depth --order 2 --fmax $fmax --dx $dx >> $logs_path/log-TB_2nd-abc_$grid_str.log;
done
