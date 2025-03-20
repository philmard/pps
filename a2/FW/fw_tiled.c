/*
 * Tiled version of the Floyd-Warshall algorithm.
 * command-line arguments: N, B
 * N = size of graph
 * B = size of tile
 * works only when N is a multiple of B
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "util.h"
#include <omp.h>

inline int min(int a, int b);
inline void FW(int **A, int K, int I, int J, int N);

int main(int argc, char **argv)
{
	int **A;
	int i,j,k;
	struct timeval t1, t2;
	double time;
	int B=64;
	int N=1024;

	if (argc != 3){
		fprintf(stdout, "Usage %s N B\n", argv[0]);
		exit(0);
	}

	N=atoi(argv[1]);
	B=atoi(argv[2]);

	A=(int **)malloc(N*sizeof(int *));
	for(i=0; i<N; i++)A[i]=(int *)malloc(N*sizeof(int));

	graph_init_random(A,-1,N,128*N);

	gettimeofday(&t1,0);

	//Iterates over each diagonal block (tile) of size BxB in the matrix
	#pragma omp parallel shared(A, N, B) private(i, j, k)
	for(k=0;k<N;k+=B){
	//	#pragma omp single{
		FW(A,k,k,k,B);
	//	}

		#pragma omp barrier		
		//Updates each tile to the left of the diagonal tile in the current row
		#pragma omp for
		for(i=0; i<k; i+=B){
			FW(A,k,i,k,B); //left
			FW(A,k,k,i,B); //above	
		}

		//Updates each tile to the right and below of the diagonal tile in the current row
		#pragma omp for
		for(i=k+B; i<N; i+=B){
			FW(A,k,i,k,B); //right
			FW(A,k,k,i,B); //below
		}

		/*
		
		//Updates each tile above the diagonal tile in the current column
		for(j=0; j<k; j+=B)
			FW(A,k,k,j,B) #pragma omp for collapse(2);
		
		//Updates each tile below the diagonal tile in the current column
		for(j=k+B; j<N; j+=B)
			FW(A,k,k,j,B);
		*/

		#pragma omp barrier
		
		/*
		//Updates each block above and to the left of the diagonal tile
		for(i=0; i<k; i+=B)
			for(j=0; j<k; j+=B)
				FW(A,k,i,j,B);

		//Updates each block above and to the right of the diagonal tile
		for(i=0; i<k; i+=B)
			for(j=k+B; j<N; j+=B)
				FW(A,k,i,j,B);

		//Updates each block below and to the left of the diagonal tile
		for(i=k+B; i<N; i+=B)
			for(j=0; j<k; j+=B)
				FW(A,k,i,j,B);

		//Updates each block below and to the right of the diagonal tile
		for(i=k+B; i<N; i+=B)
			for(j=k+B; j<N; j+=B)
				FW(A,k,i,j,B);
		*/
		
		#pragma omp for collapse(2)
    		for (i = 0; i < N; i += B) {
       			for (j = 0; j < N; j += B) {
           			if (i != k && j != k) {
               				FW(A, k, i, j, B);  // Update outer tile
           			}
        		}
    		}
		#pragma omp barrier
		
	}
	gettimeofday(&t2,0);

	time=(double)((t2.tv_sec-t1.tv_sec)*1000000+t2.tv_usec-t1.tv_usec)/1000000;
	printf("FW_TILED,%d,%d,%.4f\n", N,B,time);

/*	
	for(i=0; i<N; i++)
		for(j=0; j<N; j++) fprintf(stdout,"%d\n", A[i][j]);
	
*/	
	return 0;
}

inline int min(int a, int b)
{
	if(a<=b)return a;
	else return b;
}

inline void FW(int **A, int K, int I, int J, int N)
{
	int i,j,k;

	for(k=K; k<K+N; k++)
		for(i=I; i<I+N; i++)
			for(j=J; j<J+N; j++)
				A[i][j]=min(A[i][j], A[i][k]+A[k][j]);

}
