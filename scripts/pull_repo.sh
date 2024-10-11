git pull origin
#module load intel-oneapi-compilers-2022.0.1-gcc-7.5.0-2lzufe5
#icpc -xHost -qopenmp -O3 -I. test_TB_kernel.cpp
#make clean; make install
