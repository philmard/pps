#!/bin/bash

## Give the Job a descriptive name
#PBS -N game_of_life_performance

## Output and error files
#PBS -o game_of_life_performance.out
#PBS -e game_of_life_performance.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=1

## How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Load necessary modules
module load openmp

## Navigate to the directory with your compiled program
cd /home/parallel/parlab38/a1

## Compile the code
make
