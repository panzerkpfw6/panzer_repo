#!/bin/bash

module load intel-oneapi-compilers/2021.4.0/gcc-7.5.0-sqbobre
#module load icc/2020.2.254
module load cmake
./configure --enable-debug


#####
rm ./bin/rtm
mv -f ./CMakeCache.txt ./CMakeCache-old.txt    #Last CMakeCache.txt is saved
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install