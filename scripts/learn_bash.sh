#!/bin/bash

#for i in $(seq 1 10); do
#echo $i
#done

num_threads_max=32    #max OMP_NUM_THREADS
num_threads_min=$(expr $num_threads_max / 4)
num_threads_step=$(expr $num_threads_max / 4)
num_threads=($(seq $num_threads_min $num_threads_step $num_threads_max))

#for num in "${num_threads[@]}"; do
#  echo $num
#done

tgs_values_min=4
tgs_values_step=2
tgs_values_max=16     #controlled in the loop for tgs<num_threads
tgs_values=($(seq $tgs_values_min $tgs_values_step $tgs_values_max))

for num in "${tgs_values[@]}"; do
  echo $num
done