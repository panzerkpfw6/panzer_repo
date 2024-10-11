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
#SBATCH --job-name=test_rtm
#SBATCH --output=logs/rtm.%J.out
#SBATCH --error=logs/rtm.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee log-rtm.log) 2>&1
echo $hostname
#lscpu

###**********  workstation ***********###
###********** OPENMP PARAMETERS  ***********###
export OMP_NUM_THREADS=48
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
### export KMP_HW_SUBSET=1t
###********** MODULES & COMPILING *********###
###module load icc/2020.2.254
module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load cmake

####**********  kanary ***********###
####********** OPENMP PARAMETERS  ***********###
#export OMP_NUM_THREADS=128
#export OMP_PROC_BIND=true
#export OMP_PLACES=threads
#export OMP_NESTED='True'
#export granularity=fine
#export KMP_AFFINITY=compact
#### export KMP_HW_SUBSET=1t
####********** MODULES & COMPILING *********###
####module load icc/2020.2.254
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
#module load cmake

#####
#rm ./bin/rtm
#rm ./bin/gather
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

####*********** Delete files ************###
#rm ./data/*.raw*
#rm ./data/*.txt*
####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
timesteps=400
#nb_snap=100
nx=676;ny=676;nz=201;
first=16449;last=16450;

first=90;last=91;
dshot=2000;fmax=10;

#####*********** order 2
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1
#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1
#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last -c --fwd_steps 3 --order 2 --src_depth 5 --rcv_depth 5 --drcv 1 --dir "./data"

#####*********** order 1
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1

#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax 11 --src_depth 5 --rcv_depth 8 --drcv 1
#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax 11 --src_depth 5 --rcv_depth 8 --drcv 1
#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 --first $first --last 16454 -c --fwd_steps 3 --order 1 --src_depth 5 --rcv_depth 8 --drcv 1 --dir "./data"
## gather experiments
##./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 --first $first --last $last --src_depth 5 --drcv 1 --dir "./data" -c
##./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --dir "./data" -c
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --mode 2  --dshot 1 --first $first --last $last --src_depth 5 --drcv 1 --order 2 --fmax 8

