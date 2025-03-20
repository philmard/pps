## give the Job a descriptive name
#PBS -N small_katalipsi

## Output and error files
#PBS -o sandman_run_kmeans_critical.out
#PBS -e sandman_run_kmeans_critical.err

## How many machines should we get?
#PBS -l nodes=1:ppn=64

##How long should the job run for?
#PBS -l walltime=00:20:00

## Start
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab38/a2/a2_1


for threads in 1 2 4 8 16 32 64
do
    export OMP_NUM_THREADS=$threads
    export GOMP_CPU_AFFINITY="0-$((threads - 1))"
    echo "Running with $threads threads and CPU affinity set to $GOMP_CPU_AFFINITY"
./kmeans_omp_critical -s 32 -n 16 -c 32 -l 10
done
