#!/bin/bash

## Give the Job a descriptive name
#PBS -N KATALIPSI

## Output and error files
#PBS -o jacobi_no_conv.out
#PBS -e jacobi_no_conv.err

## Limit memory, runtime etc.
#PBS -l walltime=00:30:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=8:ppn=8

cd /home/parallel/parlab38/a4/heat_transfer/mpi

module load openmpi/1.8.3

for i in {1..3}
do
mpirun -np 64 -map-by node --mca btl tcp,self jacobi_executable 512 512 8 8
done
