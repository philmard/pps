#!/bin/bash

## Give the Job a descriptive name
#PBS -N KATALHPSH

## Output and error files
#PBS -o run_fw_sr.out
#PBS -e run_fw_sr.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=64

##How long should the job run for?
#PBS -l walltime=00:20:00

## Start 
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab38/a2/FW

nthreads=( 1 2 4 8 16 32 64 )
sizes=( 1024 2048 4096 )
block=( 256 )

for BSIZE in "${block[@]}";
do
	for size in "${sizes[@]}";
	do
		echo "--------------------------------------------------------------"
		for nthread in "${nthreads[@]}";
		do	
			echo
			export OMP_NUM_THREADS=${nthread};
			echo "Running with $nthread threads"
		##	echo "Running with $size size"
		##	echo "Running with $BSIZE bsize"
			./fw_sr ${size} ${BSIZE}
		done
	done
done
## ./fw <SIZE>
## ./fw_sr <SIZE> <BSIZE>
## ./fw_tiled <SIZE> <BSIZE>
