#!/bin/bash
## Give the Job a descriptive name
#PBS -N run_kmeans_mpi

## Output and error files
#PBS -o run_kmeans.out
#PBS -e run_kmeans.err

## Limit memory, runtime etc.
##PBS -l walltime=00:30:00

## Request 8 nodes with 8 processors (cores) per node
#PBS -l nodes=8:ppn=8

## Start
# Load the OpenMPI module (version 1.8.3)
module load openmpi/1.8.3

# Change to the directory where the kmeans executable resides
cd /home/parallel/parlab38/a4/kmeans


executable="./kmeans_mpi"
args="-c 32 -s 256 -n 16 -l 10"

for np in 1 2 4 8 16 32 64; do
    echo "--------------------------------------------------------------"
    echo "Running kmeans with ${np} processes"
##    mpirun --mca btl tcp,self -np ${np} -map-by node ${executable} ${args}
    mpirun -np $np --map-by node --mca btl self,tcp ./kmeans_mpi -s 256 -n 16 -c 32 -l 10
done

