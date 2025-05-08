#!/bin/bash
###********** OPENMP PARAMETERS ***********###
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t
###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
#module load icc/2020.2.254
module load cmake
##### COMPILATION #####
pwd


rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
#rm ./data/*ilm*
#rm ./data/*img*
#rm ./data/*sismos*
#rm ./data/*snap*
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved

##### build the project with debugging symbols
#CC=icc CXX=icpc -DCMAKE_BUILD_TYPE=Debug cmake .
#CC=icc CXX=icpc -O0 -g3 -DCMAKE_BUILD_TYPE=Debug cmake .
##### build the project
#CC=icc CXX=icpc cmake .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg" -DCMAKE_CXX_FLAGS="-pg" .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-pg -O2" -DCMAKE_CXX_FLAGS="-pg -O2" .

# debug!
#CC=icc CXX=icpc cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="-g -O0" -DCMAKE_CXX_FLAGS="-g -O0" .
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

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
tdim_arr_1st=(7 3 7)
num_wf_arr_1st=(192 20 4)

th_x_arr_1st=(2 4 4)
th_y_arr_1st=(2 2 2)
th_z_arr_1st=(1 1 1)
tdim_arr_1st=(7 3 7)
num_wf_arr_1st=(20 20 4)
export OMP_NUM_THREADS=28 #4

###*********** Experiment setup ************###
nx_arr=(  512  1024  2048  )
ny_arr=(  512  1024  2048  )
nz_arr=(  512  512   512   )
export NT_TB_1st=2001
export NT_SB_1st=2001

##### COMPILATION #####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

##### Logs directory #####
mkdir ./logs
export logs_path=./logs/test1_rtm
export logs_filename="test1_rtm_pasc.log"
rm -rf $logs_path
mkdir $logs_path

##### Run tests #####
len=${#nx_arr[@]}
#for i in $(seq 0 $len); do
for i in $(seq 1 2); do
	echo $i
	nx=${nx_arr[$i]}
	ny=${ny_arr[$i]}
	nz=${nz_arr[$i]}
	cbx=${cbx_arr[$i]}
	cby=${cby_arr[$i]}
	cbz=${cbz_arr[$i]}
	
	export logs_filename="test1_rtm_$nx.log"
	echo $logs_filename
	
#	echo "grid nx=${nx}, ny=${ny}, nz=${nz}, OMP_NUM_THREADS=${OMP_NUM_THREADS}"
#	grid_str="${nx}_${ny}_${nz}"
#	echo "cbx=${cbx}, cby=${cby}, cbz=${cbz}"
#	
#	th_x=${th_x_arr_1st[$i]}
#	th_y=${th_y_arr_1st[$i]}
#	th_z=${th_z_arr_1st[$i]}
#	t_dim=${tdim_arr_1st[$i]}
#	num_wf=${num_wf_arr_1st[$i]}
#	tgs=$((th_x * th_y*th_z))
#	echo "num_th=${OMP_NUM_THREADS},th_x=${th_x},th_y=${th_y},th_z=${th_z},num_wf=${num_wf},t_dim=${t_dim},tgs=${tgs}"
#	
#	###*********** Prepare seismograms for RTM ************###
###	srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
#	./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
#	--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
#	--tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
#	--src_depth $src_depth --order 1 --fmax $fmax --dx $dh --dy $dh --dz $dh --dt $dt --rec_sismos 1;
#
#  
#	###*********** Run RTM SB ************###
#	echo "Running RTM SB"
#	echo "Running 1st order"
#	
##	srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
#	./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st --mode 2  \
#	--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
#	--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
#	--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz >> $logs_path/$logs_filename;
#	
#	###*********** Delete RTM-related img,ilm files ************###
#	rm ./data/*img*
#	rm ./data/*ilm*
#
#	###*********** Run RTM TB ************##
#	echo "Running RTM TB"
#	echo "Running 1st order"
#	
#	#  srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 --unbuffered numactl --interleave=all \
#	./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $NT_TB_1st --tb_thread_group_size $tgs \
#	--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z \
#	--tb_t_dim $t_dim --tb_num_wf $num_wf --fwd_steps 3 --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
#	--dx $dh --dy $dh --dz $dh --dt $dt --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax $fmax >> $logs_path/$logs_filename;

done
