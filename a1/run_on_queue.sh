#!/bin/bash

## Give the Job a descriptive name
#PBS -N game_of_life_performance

## Output and error files
#PBS -o game_of_life_performance.out
#PBS -e game_of_life_performance.err

## How many machines should we get?
#PBS -l nodes=1:ppn=8

## How long should the job run for?
#PBS -l walltime=00:20:00

## Start
## Load necessary modules
module load openmp

## Navigate to the directory with your compiled program
cd /home/parallel/parlab38/a1

## Define board sizes and thread counts
board_sizes=(64 1024 4096)
thread_counts=(1 2 4 6 8)

## Loop over thread counts and board sizes
for N in "${board_sizes[@]}"; do
    for T in "${thread_counts[@]}"; do
        echo "Running Game of Life with size ${N}x${N} using ${T} threads for 1000 generations." >> performance_results.txt
        
        ## Set the number of threads
        export OMP_NUM_THREADS=$T

        ## Run the Game of Life executable and pass parameters
        ./GameOfLife $N 1000 >> performance_results.txt 2>&1  # Capture stderr as well
    done
done

echo "Performance measurements complete. Results are saved in performance_results.txt." >> performance_results.txt


