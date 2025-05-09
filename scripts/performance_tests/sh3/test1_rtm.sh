#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=workq  #
###SBATCH --mem=200G
#SBATCH --job-name=test1_rtm
#SBATCH --output=logs/test1_rtm.%J.out
#SBATCH --error=logs/test1_rtm.%J.err
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

##### Shot information #####
first=2625; # position of the source:isx=256,isy=178
last=2627; # position of the source:isx=256,isy=278
last=2625;
export src_depth=5;
export rcv_depth=8;
dshot=50;
export fmax=7;
export dh=25;
dt=0.001;

###********** Default SB, TB parameters *********###
cbx_arr=(8 16 22)
cby_arr=(6 6 46)
cbz_arr=(9999 9999 9999)

# PASC paper results
th_x_arr_1st=(16 4 4)
th_y_arr_1st=(2 2 2)
th_z_arr_1st=(1 1 1)
tdim_arr_1st=(7 3 7) # orig PASC
num_wf_arr_1st=(192 20 4)
fwd_steps=21
#fwd_steps=9

#th_x_arr_1st=(16 4 4)
#th_y_arr_1st=(2 2 4)
#th_z_arr_1st=(1 1 1)
#tdim_arr_1st=(7 3 15) # orig PASC
#num_wf_arr_1st=(192 20 4)


###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=1001
export NT_SB_1st=1001
#export NT_TB_1st=505
#export NT_SB_1st=505

#export NT_TB_1st=201
#export NT_SB_1st=201

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
export logs_path=./logs/test1_rtm
export logs_filename="test1_rtm_pasc_all_grids2.log"
#export logs_filename="test1_rtm_$nx.log"
mkdir $logs_path

rm $logs_path/*
rm $logs_path/$logs_filename;

##### Run tests #####
len=${#nx_arr[@]}
#for i in $(seq 2 $len); do
#for i in $(seq 1 2); do
for ((i=($len-1); i>=0; i--)); do
	echo $i
	nx=${nx_arr[$i]}
	ny=${ny_arr[$i]}
	nz=${nz_arr[$i]}
	cbx=${cbx_arr[$i]}
	cby=${cby_arr[$i]}
	cbz=${cbz_arr[$i]}
	cbz=9999
	
#	export logs_filename="test1_rtm_$nx.log"	
	echo $logs_filename
	
	echo "grid nx=${nx}, ny=${ny}, nz=${nz}, OMP_NUM_THREADS=${OMP_NUM_THREADS}"
	grid_str="${nx}_${ny}_${nz}"
	echo "cbx=${cbx}, cby=${cby}, cbz=${cbz}"
	
	th_x=${th_x_arr_1st[$i]}
	th_y=${th_y_arr_1st[$i]}
	th_z=${th_z_arr_1st[$i]}
	t_dim=${tdim_arr_1st[$i]}
	num_wf=${num_wf_arr_1st[$i]}
	tgs=$((th_x * th_y*th_z))
	echo "num_th=${OMP_NUM_THREADS},th_x=${th_x},th_y=${th_y},th_z=${th_z},num_wf=${num_wf},t_dim=${t_dim},tgs=${tgs}"
	
	###*********** Delete previously recorded seismograms. Clean folder. ************###
	rm -rf ./data/*sismos*
	rm ./data/*img*
	rm ./data/*ilm*
	
	###*********** Prepare seismograms for RTM ************###
	echo "Prepare seismograms for RTM"
	srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
	./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
	--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
	--tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
	--src_depth $src_depth --order 1 --fmax $fmax --dx $dh --dy $dh --dz $dh --dt $dt --rec_sismos 1 >> $logs_path/$logs_filename;
#	exit 1


	###*********** Run RTM TB ************##
	echo "Running RTM TB"
	echo "Running 1st order"
	
	srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
	./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
	--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
	--tb_t_dim $t_dim --tb_num_wf $num_wf --fwd_steps $fwd_steps --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
	--dx $dh --dy $dh --dz $dh --dt $dt --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax $fmax >> $logs_path/$logs_filename;
	###*********** Delete RTM-related img,ilm files ************###
	rm ./data/*img*
	rm ./data/*ilm*
	############################################################
  
	
	
	###*********** Run RTM SB ************###
	echo "Running RTM SB"
	echo "Running 1st order"
	
	srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
	./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st --mode 2  \
	--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
	--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
	--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/$logs_filename;
	###*********** Delete RTM-related img,ilm files ************###
	rm ./data/*img*
	rm ./data/*ilm*
#	############################################################
	
	

done
