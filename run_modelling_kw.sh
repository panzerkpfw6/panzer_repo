#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH -N 1
#SBATCH --job-name=tests_2nd_order
#SBATCH --output=logs/modelling.%J.out
#SBATCH --error=logs/modelling.%J.err
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --time=6:50:00
#SBATCH --partition=workq
#SBATCH -A k1205
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-modelling.log) 2>&1
echo $hostname
#lscpu
###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
#!!!!!!!!!
export CFLAGS=-xHost
export CXXFLAGS=-xHost

###********** MODULES & COMPILING *********###
echo preloaded_modules
module list
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5    #intel compiler: newer is better  #intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load intel/2022.1.0
#module load intel
#module load intel-oneapi/2021.4.0
#module load intel-classic/2021.4.0
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load cmake
module list

#####
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
mode=2
#timesteps=2000
timesteps=500
nx=512;ny=512;nz=512;
#nx=2048;ny=2048;nz=512;
#####***********#####***********#####***********#####***********
###*********** RUNNING SIMWAVE ************###
rm SB.log
rm TB.log
rm data/SB-final_snap.raw
rm data/TB-final_snap.raw

#numactl --interleave=all ./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --first 1301 --last 1301 --order 1 > SB1st.log
#numactl --interleave=all ./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --first 1301 --last 1301 --order 2 > SB2nd.log
./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --first 1301 --last 1301 --order 1 > SB1st.log
./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --first 1301 --last 1301 --order 2 > SB2nd.log
mv SB_lastshot_u0 data/SB-final_snap.raw

#numactl --interleave=all ./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz  --iter $timesteps  --tb_thread_group_size 4 --tb_nb_thread_groups 9 --tb_th_x 1 --tb_th_y 2 --tb_th_z 2 --tb_t_dim 3 --tb_num_wf 21 --dshot 1 --first 1301 --last 1301 -c > TB2nd.log
./bin/modeling  --verbose --n1 $nx --n2 $ny --n3 $nz  --iter $timesteps  --tb_thread_group_size 4 --tb_nb_thread_groups 9 --tb_th_x 1 --tb_th_y 2 --tb_th_z 2 --tb_t_dim 3 --tb_num_wf 21 --dshot 1 --first 1301 --last 1301 -c > TB2nd.log
mv TB_lastshot_u0 data/TB-final_snap.raw
#rm lastshot_u1

###*********** CHECKING RESULTS ***********###
./scripts/isnan data/TB-final_snap.raw
./scripts/isnan data/SB-final_snap.raw
./scripts/diff_to data/TB-final_snap.raw data/SB-final_snap.raw
#####***********#####***********#####***********#####***********



#####*********** RUNNING SIMWAVE ************###
## iterate on parameters values
#echo "CSV: size ; timesteps ; OMP_NUM_THREADS ; tgs ; ntg ; (th_x, th_y, th_z) ; t_dim ; num_wf ; mode ; tb or sb"
#size=512
## Thread configuration: SIMWAVE x=1, RACHED z=1 (innermost dimension is vectorized)
####*** temporal blocking parameters***###
##thread group : (1,2,2)  #group size: 4 #num_group : 8
#th_z=2
#th_y=2
#th_x=1
#tgs=4
#num_wf=21
#t_dim=3
#echo $(expr $OMP_NUM_THREADS / $tgs)
#echo "PROGRAM CALL: numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c"
#echo "CSV: $size ; $timesteps ; $OMP_NUM_THREADS ; $tgs ; $(expr $OMP_NUM_THREADS / $tgs) ; ($th_x, $th_y, $th_z) ; $t_dim ; $num_wf ; $mode ; tb"
####*********** TB tests ************###
## Run without First Touch
#export FIRST_TOUCH=0
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
#
### Run with First Touch
#export FIRST_TOUCH=1
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
#####*********** SB tests ************###
##numactl --interleave=all ./bin/modeling-snap_end --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
#
#echo modelling_without_snapshot_saving
#echo SB
#srun --ntasks=1 --cpus-per-task=32 --unbuffered numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
#echo SB_different
#srun --unbuffered numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##numactl --interleave=all ./bin/modeling_original --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --dshot 1 --first 1301 --last 1301
##sleep 600

