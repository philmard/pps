#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o sandman_run_kmeans.out
#PBS -e sandman_run_kmeans.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=64

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab38/a2/kmeans


for threads in 1 2 4 8 16 32 64
do
    export OMP_NUM_THREADS=$threads
    export GOMP_CPU_AFFINITY="0-$((threads - 1))"
    echo "Running with $threads threads and CPU affinity set to $GOMP_CPU_AFFINITY"
./kmeans_omp_reduction -s 256 -n 1 -c 4 -l 10
done
