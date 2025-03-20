#!/bin/bash

## Give the Job a descriptive name

#PBS -N JACOBI_MAAAAAKE

## Output and error files
#PBS -o make_jacobi.out
#PBS -e make_jacobi.err

## Limit memory, runtime etc.

#PBS -l walltime=00:05:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=1:ppn=1

cd /home/parallel/parlab38/a4/heat_transfer/mpi
module load openmpi/1.8.3
mpicc -O3 -Wall -lm mpi_jacobi.c utils.c -o jacobi_executable_no_conv
