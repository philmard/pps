#!/bin/bash
## Give the Job a descriptive name
#PBS -N make_kmeans_mpi

## Output and error files
#PBS -o make_kmeans.out
#PBS -e make_kmeans.err

## Limit memory, runtime etc.
##PBS -l walltime=00:30:00

## Request 8 nodes with 8 processors (cores) per node
#PBS -l nodes=1:ppn=8

## Start
# Load the OpenMPI module (version 1.8.3)
module load openmpi/1.8.3

cd /home/parallel/parlab38/a4/kmeans

make

