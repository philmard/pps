#!/bin/bash

## Give the Job a descriptive name
#PBS -N KATALHPSH

## How many machines should we get?
#PBS -l nodes=1:ppn=128

## How long should the job run for?
#PBS -l walltime=02:00:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab38/a2/a2_3/conc_ll

# Define the parameters for the experiments

# Number of threads
THREAD_COUNTS=(1)

# List sizes
LIST_SIZES=(1024 8192)

# Operation ratios: search-insert-delete
RATION_RATIOS=("100-0-0" "80-10-10" "20-40-40" "0-50-50")

EXECUTABLES=("x.serial")

# Iterate over all executables
for executable in "${EXECUTABLES[@]}"
do
    # Define an output file for each executable
    output_file="${executable}_output.out"

    # Iterate over all combinations of list size and operation ratio
    for list_size in "${LIST_SIZES[@]}"
    do
        for ratio in "${RATION_RATIOS[@]}"
        do
            for threads in "${THREAD_COUNTS[@]}"
            do
                # Define MT_CONF based on the number of threads
                if [ "$threads" -le 32 ]; then
                    export MT_CONF=$(seq -s, 0 $((threads - 1)))
                elif [ "$threads" -le 64 ]; then
                    export MT_CONF=$(seq -s, 0 31),$(seq -s, 0 $((threads - 33)))
                else
                    export MT_CONF=$(seq -s, 0 31),$(seq -s, 0 31),$(seq -s, 0 31),$(seq -s, 0 31)
                fi

                # Extract the operation ratios
                read -r search insert delete <<<$(echo $ratio | tr '-' ' ')

                # Log the current configuration to the output file for each executable
                echo "Running $executable with $threads threads, list size $list_size, ratio $ratio" >> $output_file

                # Ensure the arguments are passed correctly to the program
                ./$executable $list_size $search $insert $delete >> $output_file 2>&1
            done
        done
    done
done

