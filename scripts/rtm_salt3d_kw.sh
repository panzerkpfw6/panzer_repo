#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --threads-per-core=1
#SBATCH --time=24:00:00
#SBATCH --partition=7773X  # Milan-X 128
#SBATCH --job-name=test_rtm
#SBATCH --output=logs/rtm.%J.out
#SBATCH --error=logs/rtm.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading

###******** HORODATED LOG WRITING *********###
echo $hostname
#lscpu

rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
rm ./data/*ilm*
rm ./data/*img*
rm ./data/*sismos*
rm ./data/*snap*

###**********  workstation ***********###
###********** OPENMP PARAMETERS  ***********###
export OMP_NUM_THREADS=36
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
### export KMP_HW_SUBSET=1t
###********** MODULES & COMPILING *********###
###module load icc/2020.2.254
module load intel-oneapi-compilers/2022.2.1/gcc-11.3.0-k2f52ij
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
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-g -O0" .
make clean
make VERBOSE=1
make install

####*********** Delete files ************###
#rm ./data/*.raw*
#rm ./data/*.txt*
####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
#timesteps=2000
timesteps=2200
#nb_snap=100

nx=128;ny=256;nz=512;
nx=256;ny=128;nz=128;
#first=1;last=10;
first=16390;last=16506;
first=1626;last=1638;
dshot=10;
fmax=11;
#####*********** SB tests ************###
echo !!SB!!
#srun --ntasks=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --unbuffered numactl --interleave=all ./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --first 1301 --last 1301 #--nbsnap $nb_snap

#####*********** order 2
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1
#./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 2 --fmax $fmax --src_depth 5 --rcv_depth 5 --drcv 1
#./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last -c --fwd_steps 3 --order 2 --src_depth 5 --rcv_depth 5 --drcv 1 --dir "./data"

#####*********** order 1
echo "Model data for RTM"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1;

echo "Perform RTM"
./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1;

#gdb --args ./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 \
# --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax 11 --src_depth 5 --rcv_depth 8 --drcv 1;

echo "Gather images"
./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot 1 --mode 2 --first $first --last $last -c --fwd_steps 3 --order 1 --src_depth 5 --rcv_depth 8 --drcv 1 --dir "./data"


