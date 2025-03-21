#include <stdio.h>
#include <stdlib.h>

#include "kmeans.h"
#include "alloc.h"
#include "error.h"

#ifdef __CUDACC__
inline void checkCuda(cudaError_t e) {
    if (e != cudaSuccess) {
        // cudaGetErrorString() isn't always very helpful. Look up the error
        // number in the cudaError enum in driver_types.h in the CUDA includes
        // directory for a better explanation.
        error("CUDA Error %d: %s\n", e, cudaGetErrorString(e));
    }
}

inline void checkLastCudaError() {
    checkCuda(cudaGetLastError());
}
#endif

__device__ int get_tid(){
	//return 0; /* TODO: copy me from naive version... */
  return threadIdx.x + blockIdx.x * blockDim.x;
}

/* square of Euclid distance between two multi-dimensional points using column-base format */
__host__ __device__ inline static
double euclid_dist_2_transpose(int numCoords,
                    int    numObjs,
                    int    numClusters,
                    double *objects,     // [numCoords][numObjs]
                    double *clusters,    // [numCoords][numClusters]
                    int    objectId,
                    int    clusterId)
{
    int i;
    double ans=0.0;

	/* TODO: Calculate the euclid_dist of elem=objectId of objects from elem=clusterId from clusters, but for column-base format!!! */
    for(i = 0; i < numCoords; ++i) {
      double diff = (objects[objectId + i*numObjs] - clusters[clusterId + i*numClusters]);
      ans += diff * diff;
    }

    return(ans);
}

__global__ static
void find_nearest_cluster(int numCoords,
                          int numObjs,
                          int numClusters,
                          double *deviceobjects,           //  [numCoords][numObjs]
                          int *devicenewClusterSize,      // added this , [numClusters]
                          double *devicenewClusters,    // added this , [numCoords][numClusters]
                          double *deviceClusters,    //  [numCoords][numClusters]
                          int *deviceMembership,          //  [numObjs]
                          double *devdelta)
{
     extern __shared__ double shmemClusters[];

	/* TODO: copy me from shared version... */
    int t = threadIdx.x;
    int i , j;
    for(i = t; i < numClusters*numCoords; i += blockDim.x) {
        shmemClusters[i] = deviceClusters[i];
    }
    __syncthreads();

	/* Get the global ID of the thread. */
    int tid = get_tid();

    if (tid < numObjs) {

    		/* TODO: copy me from shared version... */
        int   index, i;
        double dist, min_dist;

        /* find the cluster id that has min distance to object */
        index = 0;
        /* TODO: call min_dist = euclid_dist_2(...) with correct objectId/clusterId using clusters in shmem*/
        min_dist = euclid_dist_2_transpose(numCoords , numObjs , numClusters , deviceobjects,  shmemClusters , tid , 0);

        for (i=1; i<numClusters; i++) {
            /* TODO: call dist = euclid_dist_2(...) with correct objectId/clusterId using clusters in shmem*/
            dist = euclid_dist_2_transpose(numCoords , numObjs , numClusters , deviceobjects,  shmemClusters , tid , i);
            /* no need square root */
            if (dist < min_dist) { /* find the min and its array index */
                min_dist = dist;
                index    = i;
            }
        }

        if (deviceMembership[tid] != index) {
          /* TODO: Maybe something is missing here... is this write safe? */
            //(*devdelta)+= 1.0;
            atomicAdd(devdelta , 1.0);
        }

        /* assign the deviceMembership to object objectId */
        deviceMembership[tid] = index;

    	/* TODO: additional steps for calculating new centroids in GPU? */

      atomicAdd(&devicenewClusterSize[index] , 1);

      for (j=0; j<numCoords; j++)
        atomicAdd(&devicenewClusters[j*numClusters + index] , deviceobjects[j*numObjs + tid]);
    }
}

__global__ static
void update_centroids(int numCoords,
                          int numClusters,
                          int *devicenewClusterSize,           //  [numClusters]
                          double *devicenewClusters,    //  [numCoords][numClusters]
                          double *deviceClusters)    //  [numCoords][numClusters])
{

    /* TODO: additional steps for calculating new centroids in GPU? */
    /*
    for (i=0; i<numClusters; i++) {
        for (j=0; j<numCoords; j++) {
            if (newClusterSize[i] > 0)
                dimClusters[j][i] = newClusters[j][i] / newClusterSize[i];
            newClusters[j][i] = 0.0;
        }
        newClusterSize[i] = 0;
    }
    */
    int tid = get_tid();
    int cid = tid % numClusters;
    if(tid < numCoords * numClusters) {
    	if(devicenewClusterSize[cid] > 0)
        	deviceClusters[tid] =  devicenewClusters[tid] / devicenewClusterSize[cid];

    	devicenewClusters[tid] = 0;
    }
    __syncthreads();
    devicenewClusterSize[cid] = 0;
}

//
//  ----------------------------------------
//  DATA LAYOUT
//
//  objects         [numObjs][numCoords]
//  clusters        [numClusters][numCoords]
//  dimObjects      [numCoords][numObjs]
//  dimClusters     [numCoords][numClusters]
//  newClusters     [numCoords][numClusters]
//  deviceObjects   [numCoords][numObjs]
//  deviceClusters  [numCoords][numClusters]
//  ----------------------------------------
//
/* return an array of cluster centers of size [numClusters][numCoords]       */
void kmeans_gpu(	double *objects,      /* in: [numObjs][numCoords] */
		               	int     numCoords,    /* no. features */
		               	int     numObjs,      /* no. objects */
		               	int     numClusters,  /* no. clusters */
		               	double   threshold,    /* % objects change membership */
		               	long    loop_threshold,   /* maximum number of iterations */
		               	int    *membership,   /* out: [numObjs] */
						double * clusters,   /* out: [numClusters][numCoords] */
						int blockSize)
{
    double timing = wtime(), timing_internal, timer_min = 1e42, timer_max = 0;
    double timer1 = 0 , timer2 = 0 , temp;
	int    loop_iterations = 0;
    int      i, j, index, loop=0;
    double  delta = 0, *dev_delta_ptr;          /* % of objects change their clusters */
    /* TODO: Copy me from transpose version*/
    double  **dimObjects = (double**) calloc_2d(numCoords , numObjs , sizeof(double)); //calloc_2d(...) -> [numCoords][numObjs]
    double  **dimClusters = (double**) calloc_2d(numCoords , numClusters , sizeof(double)); //calloc_2d(...) -> [numCoords][numClusters]
    double  **newClusters = (double**) calloc_2d(numCoords , numClusters , sizeof(double));  //calloc_2d(...) -> [numCoords][numClusters]

    printf("\n|-----------Full-offload GPU Kmeans------------|\n\n");

    /* TODO: Copy me from transpose version*/
    for(i = 0; i < numCoords; ++i) {
                for(j = 0; j < numObjs; ++j) {
                        dimObjects[i][j] = objects[i + j*numCoords];
                }
    }

    double *deviceObjects;
    double *deviceClusters, *devicenewClusters;
    int *deviceMembership;
    int *devicenewClusterSize; /* [numClusters]: no. objects assigned in each new cluster */

    /* pick first numClusters elements of objects[] as initial cluster centers*/
    for (i = 0; i < numCoords; i++) {
        for (j = 0; j < numClusters; j++) {
            dimClusters[i][j] = dimObjects[i][j];
        }
    }

    /* initialize membership[] */
    for (i=0; i<numObjs; i++) membership[i] = -1;

    timing = wtime() - timing;
    printf("t_alloc: %lf ms\n\n", 1000*timing);
    timing = wtime();
    const unsigned int numThreadsPerClusterBlock = (numObjs > blockSize)? blockSize: numObjs;
    const unsigned int numClusterBlocks = ((numObjs-1) / numThreadsPerClusterBlock) + 1; /* TODO: Calculate Grid size, e.g. number of blocks. */
	/*	Define the shared memory needed per block.
    	- BEWARE: We can overrun our shared memory here if there are too many
    	clusters or too many coordinates!
    	- This can lead to occupancy problems or even inability to run.
    	- Your exercise implementation is not requested to account for that (e.g. always assume deviceClusters fit in shmemClusters */
    const unsigned int clusterBlockSharedDataSize = numClusters * numCoords * sizeof(double);

    cudaDeviceProp deviceProp;
    int deviceNum;
    cudaGetDevice(&deviceNum);
    cudaGetDeviceProperties(&deviceProp, deviceNum);

    if (clusterBlockSharedDataSize > deviceProp.sharedMemPerBlock) {
        error("Your CUDA hardware has insufficient block shared memory to hold all cluster centroids\n");
    }

    checkCuda(cudaMalloc(&deviceObjects, numObjs*numCoords*sizeof(double)));
    checkCuda(cudaMalloc(&deviceClusters, numClusters*numCoords*sizeof(double)));
    checkCuda(cudaMalloc(&devicenewClusters, numClusters*numCoords*sizeof(double)));
    checkCuda(cudaMalloc(&devicenewClusterSize, numClusters*sizeof(int)));
    checkCuda(cudaMalloc(&deviceMembership, numObjs*sizeof(int)));
    checkCuda(cudaMalloc(&dev_delta_ptr, sizeof(double)));

    timing = wtime() - timing;
    printf("t_alloc_gpu: %lf ms\n\n", 1000*timing);
    timing = wtime();

    checkCuda(cudaMemcpy(deviceObjects, dimObjects[0],
              numObjs*numCoords*sizeof(double), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(deviceMembership, membership,
              numObjs*sizeof(int), cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(deviceClusters, dimClusters[0],
                  numClusters*numCoords*sizeof(double), cudaMemcpyHostToDevice));
    checkCuda(cudaMemset(devicenewClusterSize, 0, numClusters*sizeof(int)));

    //added this
    checkCuda(cudaMemset(devicenewClusters , 0.0 , numClusters*numCoords*sizeof(double)));



    free(dimObjects[0]);

    timing = wtime() - timing;
    printf("t_get_gpu: %lf ms\n\n", 1000*timing);
    timing = wtime();

    do {
        timing_internal = wtime();
        checkCuda(cudaMemset(dev_delta_ptr, 0, sizeof(double)));
	
	//temp = wtime();
		//printf("Launching find_nearest_cluster Kernel with grid_size = %d, block_size = %d, shared_mem = %d KB\n", numClusterBlocks, numThreadsPerClusterBlock, clusterBlockSharedDataSize/1000);
        /* TODO: change invocation if extra parameters needed */
        find_nearest_cluster
            <<< numClusterBlocks, numThreadsPerClusterBlock, clusterBlockSharedDataSize >>>
            (numCoords, numObjs, numClusters,
             deviceObjects, devicenewClusterSize, devicenewClusters, deviceClusters, deviceMembership, dev_delta_ptr);

        cudaDeviceSynchronize(); checkLastCudaError();
		//printf("Kernels complete for itter %d, updating data in CPU\n", loop);
	
	//timer1 += wtime() - temp;

    	/* TODO: Copy dev_delta_ptr to &delta */
        checkCuda(cudaMemcpy(&delta , dev_delta_ptr , sizeof(double) , cudaMemcpyDeviceToHost));


      /* TODO: calculate dim for "update_centroids" and fire it  */
      // we will have numCoords * numClusters threads so that each thread will be a one-dimensional coordinate of a kernel
     	const unsigned int update_centroids_block_sz = (numCoords* numClusters > blockSize) ? blockSize: numCoords* numClusters;  /* TODO: can use different blocksize here if deemed better */
     	const unsigned int update_centroids_dim_sz =  (((numCoords*numClusters)-1) / update_centroids_block_sz) + 1;
	
	//temp = wtime();
     	update_centroids<<< update_centroids_dim_sz, update_centroids_block_sz, 0 >>>
            (numCoords, numClusters, devicenewClusterSize, devicenewClusters, deviceClusters);
        cudaDeviceSynchronize(); checkLastCudaError();
	//timer2 += wtime() - temp;

        delta /= numObjs;
       	//printf("delta is %f - ", delta);
        loop++;
        //printf("completed loop %d\n", loop);
		timing_internal = wtime() - timing_internal;
		if ( timing_internal < timer_min) timer_min = timing_internal;
		if ( timing_internal > timer_max) timer_max = timing_internal;
	} while (delta > threshold && loop < loop_threshold);

    checkCuda(cudaMemcpy(membership, deviceMembership,
                 numObjs*sizeof(int), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(dimClusters[0], deviceClusters,
                 numClusters*numCoords*sizeof(double), cudaMemcpyDeviceToHost));

	for (i=0; i<numClusters; i++) {
		for (j=0; j<numCoords; j++) {
		    clusters[i*numCoords + j] = dimClusters[j][i];
		}
	}

    timing = wtime() - timing;
    printf("nloops = %d  : total = %lf ms\n\t-> t_loop_avg = %lf ms\n\t-> t_loop_min = %lf ms\n\t-> t_loop_max = %lf ms\n\n|-------------------------------------------|\n",
    	loop, 1000*timing, 1000*timing/loop, 1000*timer_min, 1000*timer_max);
	
	//printf("First kernel: %lf \n", 1000*timer1);
	//printf("Second kernel: %lf \n", 1000*timer2);

	char outfile_name[1024] = {0};
	sprintf(outfile_name, "Execution_logs/silver1-V100_Sz-%lu_Coo-%d_Cl-%d.csv", numObjs*numCoords*sizeof(double)/(1024*1024), numCoords, numClusters);
	FILE* fp = fopen(outfile_name, "a+");
	if(!fp) error("Filename %s did not open succesfully, no logging performed\n", outfile_name);
	fprintf(fp, "%s,%d,%lf,%lf,%lf\n", "All_GPU", blockSize, timing/loop, timer_min, timer_max);
	fclose(fp);

    checkCuda(cudaFree(deviceObjects));
    checkCuda(cudaFree(deviceClusters));
    checkCuda(cudaFree(devicenewClusters));
    checkCuda(cudaFree(devicenewClusterSize));
    checkCuda(cudaFree(deviceMembership));

    return;
}

