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
#SBATCH --job-name=test_rtm_2nd
#SBATCH --output=logs/rtm.%J.out
#SBATCH --error=logs/rtm.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-modeling_.log) 2>&1
echo $hostname
#lscpu
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=48
#export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
#export KMP_HW_SUBSET=1t

###********** MODULES *********###
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load icc/2020.2.254
module load cmake

##### COMPILATION #####
#cd ..
rm ./bin/modeling
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
#CC=icc CXX=icpc -DCMAKE_BUILD_TYPE=Debug cmake .
CC=icc CXX=icpc cmake .

make clean
make VERBOSE=1
make install

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
mode=2
export TIME_SB_1st=504
#export TIME_SB_1st=2
export TIME_TB_1st=505
#export TIME_TB_1st=1000
#export TIME_TB_1st=512
#nx=512;ny=512;nz=512;
#nx=2048;ny=2048;nz=512;
nx=128;ny=256;nz=512;
####*********** Directories ************###
fld=./logs2
mkdir $fld
rm $fld/*

#####*********** SB tests ************########***********
## choice of shot: the shot index corresponds to the receiver index:
## s->shots[idx]->srcidx = s->rcv[ir] + s->drcv / 2;
## x,y indexing in acquisition.txt file is relative to the boundaries of computation domain().
export shot=16447
rm snapshot_504

./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode $mode --dshot 1 --first $shot --last $shot --src_depth 256 --drcv 1 --order 1 --fmax 8
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_SB_1st --mode 2  --dshot 1 --first $shot --last $shot --src_depth 256 --drcv 1 --order 2 --fmax 8
#./bin/modeling --verbose --n1 128 --n2 256 --n3 512 --iter $TIME_SB_1st --mode 2  --dshot 1 --first $shot --last $shot --src_depth 256 --drcv 1 --order 1
#####*********** TB tests ************########***********
echo !!TB with parameters!!
#export FIRST_TOUCH=1m
#srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
th_x=8
th_y=2
th_z=1
tgs=16
t_dim=7
num_wf=64

echo "TB_mod_1st"
#./bin/modeling --verbose --n1 $nx  --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2  --dshot 1 --first $shot --last $shot  --fwd_steps 3 -c --order 1 --fmax 8
echo "TB_mod_2nd"
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2  --dshot 1 --first $shot --last $shot --fwd_steps 3 -c --order 2 --fmax 8
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $TIME_TB_1st --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf --mode 2  --dshot 1 --first $shot --last $shot --fwd_steps 3 -c
#./bin/modeling --verbose --n1 128 --n2 256 --n3 512 --iter 505 --tb_thread_group_size 16 --tb_nb_thread_groups 1 --tb_th_x 8 --tb_th_y 2 --tb_th_z 1 --tb_t_dim 7 --tb_num_wf 64 --mode 2 --dshot 1 --first $shot --last $shot -c
#./bin/modeling --verbose --n1 128 --n2 256 --n3 512 --iter 505 --tb_thread_group_size 1 --tb_nb_thread_groups 1 --tb_th_x 1 --tb_th_y 1 --tb_th_z 1 --tb_t_dim 7 --tb_num_wf 64 --mode 2 --dshot 1 --first $shot --last $shot -c

### order 1
#./bin/modeling --verbose --n1 128 --n2 256 --n3 512 --iter $TIME_TB_1st --tb_thread_group_size 16 --tb_nb_thread_groups 3 --tb_th_x 8 --tb_th_y 2 --tb_th_z 1 --tb_t_dim 7 --tb_num_wf 64 --mode 2 --dshot 1 --first 16447 --last 16447 -c --fwd_steps 3 --order 1 --fmax 8 --src_depth 256 --drcv 1
### order 2
#./bin/modeling --verbose --n1 128 --n2 256 --n3 512 --iter $TIME_TB_1st --tb_thread_group_size 16 --tb_nb_thread_groups 3 --tb_th_x 8 --tb_th_y 2 --tb_th_z 1 --tb_t_dim 7 --tb_num_wf 64 --mode 2 --dshot 1 --first 16447 --last 16447 -c --fwd_steps 3 --order 2 --fmax 8 --src_depth 256 --drcv 1

