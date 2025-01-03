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
#SBATCH --job-name=reproduce_sb
#SBATCH --output=./logs/reproduce_sb.%J.out
#SBATCH --error=./logs/reproduce_sb.%J.err
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread    # don't use hyperthreading
########################################
###******** HORODATED LOG WRITING *********###
exec > >(while read line; do echo "$(date): $line"; done | tee test2.log) 2>&1
echo $hostname
#lscpu

###********** OPENMP PARAMETERS ***********###
#export OMP_NUM_THREADS=48
export OMP_NUM_THREADS=128
export OMP_PROC_BIND=true
export OMP_PLACES=threads
export OMP_NESTED='True'
export granularity=fine
export KMP_AFFINITY=compact
export KMP_HW_SUBSET=1t

###********** MODULES *********###
#module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
module load icc/2020.2.254
module load cmake

####################################################
export TIME_TB_2nd=530 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_SB_2nd=529 #@pavel in SB the nt should one time less than in correponding TB.
export TIME_TB_1st=537 #@pavel in TB source injection starts from second time sample (Nothing happens for one dt).This is code feature.
export TIME_SB_1st=536 #@pavel in SB the nt should one time less than in correponding TB.
pwd

##################### COMPILE
rm SB_1st-abc.out
rm SB_2nd_abc.out;
rm TB_1st_abc.out;
rm TB_2nd_abc.out;
cd data


cd ../SB_abc;
x=10;y=22;
sed -i "s/const Myint _CB_SIZE_X = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_X = $x  ;/g" SB_kernel.h;
sed -i "s/const Myint _CB_SIZE_Y = [[:digit:]]\+[ ]*;/const Myint _CB_SIZE_Y = $y;/g" SB_kernel.h;
icpc -xHost -qopenmp -march=core-avx2 -mtune=core-avx2 -O3 -I. test_SB_kernel.cpp -o ../SB_1st-abc.out
cd ..
#######################################
export logs_path="./logs/reproduce_sb"
mkdir "$logs_path"

nx=512; ny=512; nz=512;
grid_str="${nx}_${ny}_${nz}_${x}_${y}"
################  1st   ################
srun --nodes=1 --cpus-per-task=$OMP_NUM_THREADS --hint=nomultithread --threads-per-core=1 ./SB_1st-abc.out $nx $ny $nz $TIME_SB_1st 3000 0 0 0 0 >> $logs_path/log-SB_1st-abc_$grid_str.log
##################################################

