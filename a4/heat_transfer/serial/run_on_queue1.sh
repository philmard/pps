ve the Job a descriptive name
#PBS -N jacobi

## Output and error files
#PBS -o jacobi_no_conv_1_proc.out
#PBS -e jacobi_no_conv_1_proc.err

## Limit memory, runtime etc.
#PBS -l walltime=00:45:00

## How many nodes:processors_per_node should we get?
#PBS -l nodes=8:ppn=8

sizes='2048 4096 6144'
##sizes='2048'
procs=(1)
procs_dimy=(1)
procs_dimx=(1)

cd /home/parallel/parlab38/a4/mpi

module load openmpi/1.8.3

for size in $sizes
do
    for j in ${!procs[@]}
    do 
        mpirun -np ${procs[$j]} -map-by node --mca btl tcp,self jacobi_test $size $size ${procs_dimx[$j]} ${procs_dimy[$j]}
    done
done
