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
#SBATCH --output=logs/rtm_2layer_.%J.out
#SBATCH --error=logs/rtm_2layer_.%J.err
#SBATCH --cpus-per-task=192
#SBATCH --hint=nomultithread    # don't use hyperthreading

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

##### CREATE FOLDERS #####
mkdir data

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
##########CC=icc CXX=icpc cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="-g -O0" -DCMAKE_CXX_FLAGS="-g -O0" .
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install
export logs_path=./logs/rtm
mkdir $logs_path

################################
echo "Test 1. compare_rtm_result_for_SB_TB"
#@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_TB_1st=4000 
export TIME_SB_1st=5500
dt=0.001;fmax=7;
### grid size 256*256*256
nx=256;ny=256;nz=256;dh=25;
nx=200;ny=256;nz=160;dh=25;

x=2; y=2; z=1; t=7; w=20; tgs=4;

#x=1; y=1; z=1; t=7; w=20; tgs=4;

cbx=64;cby=22;cbz=9999;

first=513; # position of the source:isx=100,isy=50
#last=513; # position of the source:isx=160,isy=140
last=517; # position of the source:isx=100,isy=250
#first=513; last=514;
dshot=50;
export src_depth=5;
export rcv_depth=8;
##################### SB RTM workflow	#####################
export OMP_NUM_THREADS=192 #4
mkdir ./data/rtm_sb
rm ./data/rtm_sb/*

srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz  >> $logs_path/log_model_2layer.log;
#
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot \
--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz  >> $logs_path/log_rtm_2layer.log;


### move data to SB folder
mv ./data/img* ./data/rtm_sb/
mv ./data/ilm* ./data/rtm_sb/
mv ./data/sismos* ./data/rtm_sb/
###
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 \
./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot --dir "./data/rtm_sb"  \
--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz -c >> $logs_path/log_gather_2layer.log;

exit 0

##################### TB RTM workflow	##################### 
export OMP_NUM_THREADS=4 #4
#mkdir ./data/rtm_tb
#rm ./data/rtm_tb/*
# gdb --args 
#
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
#--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
#--tb_t_dim $t --tb_num_wf $w --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
#--dx $dh --dy $dh --dz $dh --dt $dt --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax $fmax;

#gdb --batch --ex "run" --ex "bt" --args ./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
#gdb --args ./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \


#gdb --args ./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
./bin/rtm --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs \
--tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $x --tb_th_y $y --tb_th_z $z \
--tb_t_dim $t --tb_num_wf $w --fwd_steps 3 --mode 2 --drcv 1 --dshot $dshot --first $first --last $last -c \
--dx $dh --dy $dh --dz $dh --dt $dt --src_depth $src_depth --rcv_depth $rcv_depth --order 1 --fmax $fmax;
 
### move data to TB folder
#mv ./data/img* ./data/rtm_tb/
#mv ./data/ilm* ./data/rtm_tb/
#mv ./data/sismos* ./data/rtm_tb/
###

#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  \
#--first $first --last $last --src_depth $src_depth --rcv_depth $rcv_depth \
#--dx $dh --dy $dh --dz $dh --dt $dt --drcv 1 --dshot $dshot --dir "./data/rtm_tb"  \
#--order 1 --fmax $fmax --cbx $cbx --cby $cby --cbz $cbz -c;

######################################################### compare data
#./scripts_useful/diff_to ./snapshot_TB1st_505 ./snapshot_SB1st_505
#./scripts_useful/diff_to ./data/rtm_sb/sismos_${shot}.raw ./data/rtm_tb/sismos_${shot}.raw
#
