#!/bin/bash
# ... (previous parts unchanged up to COMPILATION)

##### COMPILATION #####
export CFLAGS="-O3 -qopenmp -g"
export CXXFLAGS="-O3 -qopenmp -g"
export FFLAGS="-O3 -qopenmp -g"
CC=icc CXX=icpc cmake .
make clean
make VERBOSE=1
make install

#####################
echo "Running SB"
echo "Running 1st order"
./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx >> $logs_path/log-SB_1st-abc_2_$grid_str.log

#####################
echo "Profiling with perf..."
sudo sysctl -w kernel.perf_event_paranoid=1
sudo perf record -F 999 -g ./bin/modeling --verbose --n1 $nx --n2 $ny --n3 $nz --iter $NT_SB_1st \
--mode 2 --drcv 1 --dshot 1 --first $shot --last $shot --src_depth $src_depth --order 1 --fmax $fmax \
--dx $dx
sudo chown plotnips:plotnips perf.data  # Fix ownership
perf report > $logs_path/perf_report_$grid_str.txt