#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=workq  #
#SBATCH --job-name=test1_modeling
#SBATCH --output=logs/test1_modeling.%J.out
#SBATCH --error=logs/test1_modeling.%J.err
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


###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#module load icc/2020.2.254
module load cmake
##### COMPILATION #####
pwd
#rm *snapshot*
#exit 0
rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved

##### build the project with debugging symbols
#CC=icc CXX=icpc -DCMAKE_BUILD_TYPE=Debug cmake .
#CC=icc CXX=icpc -O0 -g3 -DCMAKE_BUILD_TYPE=Debug cmake .

##### build the project
#CC=icc CXX=icpc cmake .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg" -DCMAKE_CXX_FLAGS="-pg" .
CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg -O2" -DCMAKE_CXX_FLAGS="-pg -O2" .
make clean
make VERBOSE=1
make install


##### Shot information #####
export shot=32896;  # position of the source in x,y coordinates.check ./data/acquisition.txt
export src_depth=256;
export fmax=8;
export dh=10;
export dt=0.001;

###********** Default SB, TB parameters *********###
cbx_arr=(8  16 22)
cby_arr=(6  6  46)
cbz_arr=(9999  9999  9999)

## My recent parameter search
th_x_arr_1st=(16 16 16)
th_y_arr_1st=(2 2 2)
th_z_arr_1st=(1 1 1)
tdim_arr_1st=(7 7 7)
num_wf_arr_1st=(192 192 192)

## PASC paper results
#th_x_arr_1st=(16 4 4)
#th_y_arr_1st=(2 2 2)
#th_z_arr_1st=(1 1 1)
#tdim_arr_1st=(7 3 7)
#num_wf_arr_1st=(192 20 4)

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=505
export NT_SB_1st=505
export NT_TB_2nd=502
export NT_SB_2nd=505


##### Logs directory #####
mkdir ./logs
export logs_path=./logs/test1_modeling
######rm -rf $logs_path
mkdir $logs_path

##### Run tests #####
len=${#nx_arr[@]}
#for i in $(seq 0 $len); do
for i in $(seq 0 2); do
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
#  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
  gdb --args ./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_SB_1st \
  --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
  --dx $dh --dy $dh --dz $dh --dt $dt  \
  --cbx $cbx --cby $cby --cbz $cbz ;
  
#  >> $logs_path/test1_fwdSB.log;

  ###*********** TB ************##
#  echo "Running TB"
#  echo "Running 1st order"
#  th_x=${th_x_arr_1st[$i]}
#  th_y=${th_y_arr_1st[$i]}
#  th_z=${th_z_arr_1st[$i]}
#  t_dim=${tdim_arr_1st[$i]}
#  num_wf=${num_wf_arr_1st[$i]}
#  tgs=$((th_x * th_y*th_z))
#  echo "num_th=${OMP_NUM_THREADS},th_x=${th_x},th_y=${th_y},th_z=${th_z},num_wf=${num_wf},t_dim=${t_dim},tgs=${tgs}"
#
#  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
#  ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
#  --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
#  --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot 1 --first $shot --last $shot -c \
#  --src_depth $src_depth --order 1 --fmax $fmax --dx $dx >> $logs_path/test1_fwd2.log;
done
