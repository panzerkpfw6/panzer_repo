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
lscpu

rm ./bin/modeling
rm ./bin/rtm
rm ./bin/gather
#rm ./data/*ilm*
#rm ./data/*img*
#rm ./data/*sismos*
#rm ./data/*snap*

###**********  workstation ***********###
###********** OPENMP PARAMETERS  ***********###
export OMP_NUM_THREADS=44
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

####********** MODULES & COMPILING *********###
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
#CC=icc CXX=icpc cmake -DCMAKE_C_FLAGS="-g -O0" .
make clean
make VERBOSE=1
make install

####*********** RUNNING RTM ************###
###********** mode, grid, time steps ***********###
timesteps=4000;	dt=0.001;
nx=676;ny=676;nz=201;dh=25;
##### Profile x=310. Salt3D. no mistake
#first=20957;last=21023;
#### Profile x=?. Salt3D. no mistake
first=1;last=100;
first=61;last=61;
dshot=4568;
fmax=11;
cbx=64;cby=22;cbz=9999;

#####*********** order 1
#echo "Model data for RTM"
#first=20468;last=20507;
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
# --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1;
#exit 1;

#####*********** Back up data 
### Configuration
#DATA_DIR="./data"
#ORIG_DATA_DIR="./data/orig_data"
#nproc=4
#mkdir $DATA_DIR
#mkdir $ORIG_DATA_DIR
#PATTERN="*sismos*"
## Ensure destination directory exists
#mkdir -p "$ORIG_DATA_DIR"
## Find all files matching the pattern and copy them in parallel
#find "$DATA_DIR" -maxdepth 1 -type f -name "$PATTERN" | parallel -j "$(nproc)" cp {} "$ORIG_DATA_DIR"
#echo "Backup complete."
########exit 0;
### remove_direct_wave from the modeled real data
#scp ./data/*sismos* ./data/orig_data
##python ./scripts_useful/python/remove_direct_wave.py
########################
echo "Model data for RTM. salt3d."
#./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 --first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 --cbx $cbx --cby $cby --cbz $cbz;

./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
--dx $dh --dy $dh --dz $dh --dt $dt  \
--cbx $cbx --cby $cby --cbz $cbz;

echo "Perform RTM"
./bin/rtm --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
--first $first --last $last --fwd_steps 3 --order 1 --fmax $fmax --src_depth 5 --rcv_depth 8 --drcv 1 \
--dx $dh --dy $dh --dz $dh \
--cbx $cbx --cby $cby --cbz $cbz;


echo "Gather images"
./bin/gather --verbose --n1 $nx --n2 $ny --n3 $nz --iter $timesteps --dshot $dshot --mode 2 \
 --first $first --last $last -c --fwd_steps 3 --order 1 --src_depth 5 --rcv_depth 8 --drcv 1 --dir "./data" \
--cbx $cbx --cby $cby --cbz $cbz;

