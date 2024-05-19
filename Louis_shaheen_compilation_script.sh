#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH -N 1
#SBATCH --job-name=test_simwave_louis
#SBATCH --output=logs/modelling.%J.out
#SBATCH --error=logs/modelling.%J.err
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --time=3:00:00
#SBATCH --partition=workq
#SBATCH -A k1205
#SBATCH --hint=nomultithread    # don't use hyperthreading
###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-simwave_perf.log) 2>&1

###*** PARAMETERS TO ITERATE A SWEEP ON ***###
# Values to iter for num_threads | tgs | size |
num_threads_max=32    #max OMP_NUM_THREADS
num_threads_min=$(expr $num_threads_max / 4)
num_threads_step=$(expr $num_threads_max / 4)
num_threads=($(seq $num_threads_min $num_threads_step $num_threads_max))

tgs_values_min=4
tgs_values_step=2
tgs_values_max=16     #controlled in the loop for tgs<num_threads
tgs_values=($(seq $tgs_values_min $tgs_values_step $tgs_values_max))

size_values=(512 1024 2048)
timesteps=1000
num_wf=21
t_dim=3
mode=2

###********** OPENMP PARAMETERS ***********###
export OMP_NUM_THREADS=$num_threads_max
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'

###********** MODULES & COMPILING *********###
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5    #intel compiler: newer is better  #intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load intel-oneapi/2021.4.0
module load cmake
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make
make install

###*********** RUNNING SIMWAVE ************###
# iterate on parameters values
echo "CSV: size ; timesteps ; OMP_NUM_THREADS ; tgs ; ntg ; (th_x, th_y, th_z) ; t_dim ; num_wf ; mode ; tb or sb"