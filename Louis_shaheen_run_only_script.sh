#!/bin/bash
###********* SLURM CONFIGURATION **********###
#SBATCH -N 1
#SBATCH --job-name=test_simwave_louis
#SBATCH --output=logs/modelling.%J.out
#SBATCH --error=logs/modelling.%J.err
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=pavel.plotnitskii@kaust.edu.sa
#SBATCH --time=24:00:00
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
module load intel-oneapi/2021.4.0
make
###*********** RUNNING SIMWAVE ************###
# iterate on parameters values
echo "CSV: size ; timesteps ; OMP_NUM_THREADS ; tgs ; ntg ; (th_x, th_y, th_z) ; t_dim ; num_wf ; mode ; tb or sb"
for num in "${num_threads[@]}"; do
  export OMP_NUM_THREADS=$num
  for tgs in "${tgs_values[@]}"; do
    # Skip if thread_group_size is greater than OMP_NUM_THREADS
    if (( tgs > num )); then
      continue
    fi
    for size in "${size_values[@]}"; do
      # Thread configuration: SIMWAVE x=1, RACHED z=1 (innermost dimension is vectorized)
      th_z=$(expr $tgs / 2)
      th_y=2
      th_x=1

      echo "PROGRAM CALL: numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c"
      echo "CSV: $size ; $timesteps ; $OMP_NUM_THREADS ; $tgs ; $(expr $OMP_NUM_THREADS / $tgs) ; ($th_x, $th_y, $th_z) ; $t_dim ; $num_wf ; $mode ; tb"
      # Run without First Touch
      export FIRST_TOUCH=0
      numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c

      # Run with First Touch
      export FIRST_TOUCH=1
      numactl --interleave=all ./bin/modeling --verbose --n1 $size  --n2 $size --n3 512 --iter $timesteps --tb_thread_group_size $tgs --tb_nb_thread_groups $(expr $OMP_NUM_THREADS / $tgs) --tb_th_x $th_x --tb_th_y $th_y --tb_th_z $th_z --tb_t_dim $t_dim --tb_num_wf $num_wf   --mode $mode  --dshot 1 --first 1301 --last 1301  --fwd_steps 3 -c
    done
  done
done