#!/bin/bash

## Give the Job a descriptive name
#PBS -N KATALHPSH

## Output and error files
#PBS -o run_kmeans_general.out
#PBS -e run_kmeans_general.err

## How many machines should we get?
#PBS -l nodes=1:ppn=64

## How long should the job run for?
#PBS -l walltime=00:30:00

## Start
## Load OpenMP module and navigate to the directory
module load openmp
cd /home/parallel/parlab38/a2/a2_1

## List of lock mechanisms
lock_mechanisms=(
    "kmeans_omp_nosync_lock"
    "kmeans_omp_pthread_mutex_lock"
    "kmeans_omp_pthread_spin_lock"
    "kmeans_omp_tas_lock"
    "kmeans_omp_ttas_lock"
    "kmeans_omp_array_lock"
    "kmeans_omp_clh_lock"
)

## Loop through lock mechanisms and threads
for lock in "${lock_mechanisms[@]}"
do
    # Create unique output and error files for each mechanism
    out_file="run_${lock}.out"
    err_file="run_${lock}.err"

    for threads in 1 2 4 8 16 32 64
    do
        export OMP_NUM_THREADS=$threads
        export GOMP_CPU_AFFINITY="0-$((threads - 1))"
        echo "Running $lock with $threads threads and CPU affinity set to $GOMP_CPU_AFFINITY" >> $out_file
        ./$lock -s 32 -n 16 -c 32 -l 10 >> $out_file 2>> $err_file
    done
done

