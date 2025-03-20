#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "utils.h"

int main(int argc, char ** argv) {
    int rank,size;
    int global[2],local[2]; //global matrix dimensions and local matrix dimensions (2D-domain, 2D-subdomain)
    int global_padded[2];   //padded global matrix dimensions (if padding is not needed, global_padded=global)
    int grid[2];            //processor grid dimensions
    int i,j,t;
    int global_converged=0,converged=0; //flags for convergence, global and per process
    MPI_Datatype dummy;     //dummy datatype used to align user-defined datatypes in memory
    double omega; 			//relaxation factor - useless for Jacobi

    struct timeval tts,ttf,tcs,tcf,tes,tef,tss,tsf;   //Timers: total-> tts,ttf, computation -> tcs,tcf
    double ttotal=0,tcomp=0,tcomm=0,tconv=0,total_time,comp_time,comm_time,conv_time;
    
    double ** U, ** u_current, ** u_previous, ** swap; //Global matrix, local current and previous matrices, pointer to swap between current and previous
    

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    //----Read 2D-domain dimensions and process grid dimensions from stdin----//

    if (argc!=5) {
        fprintf(stderr,"Usage: mpirun .... ./exec X Y Px Py");
        exit(-1);
    }
    else {
        global[0]=atoi(argv[1]);
        global[1]=atoi(argv[2]);
        grid[0]=atoi(argv[3]);
        grid[1]=atoi(argv[4]);
    }

    //----Create 2D-cartesian communicator----//
	//----Usage of the cartesian communicator is optional----//

    MPI_Comm CART_COMM;         //CART_COMM: the new 2D-cartesian communicator
    int periods[2]={0,0};       //periods={0,0}: the 2D-grid is non-periodic
    int rank_grid[2];           //rank_grid: the position of each process on the new communicator
		
    MPI_Cart_create(MPI_COMM_WORLD,2,grid,periods,0,&CART_COMM);    //communicator creation
    MPI_Cart_coords(CART_COMM,rank,2,rank_grid);	                //rank mapping on the new communicator

    //----Compute local 2D-subdomain dimensions----//
    //----Test if the 2D-domain can be equally distributed to all processes----//
    //----If not, pad 2D-domain----//
    
    for (i=0;i<2;i++) {
        if (global[i]%grid[i]==0) {
            local[i]=global[i]/grid[i];
            global_padded[i]=global[i];
        }
        else {
            local[i]=(global[i]/grid[i])+1;
            global_padded[i]=local[i]*grid[i];
        }
    }

	//Initialization of omega
    omega=2.0/(1+sin(3.14/global[0]));

    //----Allocate global 2D-domain and initialize boundary values----//
    //----Rank 0 holds the global 2D-domain----//
    if (rank==0) {
        U=allocate2d(global_padded[0],global_padded[1]);   
        init2d(U,global[0],global[1]);
    }

    //----Allocate local 2D-subdomains u_current, u_previous----//
    //----Add a row/column on each size for ghost cells----//

    u_previous=allocate2d(local[0]+2,local[1]+2);
    u_current=allocate2d(local[0]+2,local[1]+2);   
    
    //----Distribute global 2D-domain from rank 0 to all processes----//
         
 	//----Appropriate datatypes are defined here----//
	/*****The usage of datatypes is optional*****/
    
    //----Datatype definition for the 2D-subdomain on the global matrix----//

    MPI_Datatype global_block;
    MPI_Type_vector(local[0],local[1],global_padded[1],MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&global_block);
    MPI_Type_commit(&global_block);

    //----Datatype definition for the 2D-subdomain on the local matrix----//

    MPI_Datatype local_block;
    MPI_Type_vector(local[0],local[1],local[1]+2,MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&local_block);
    MPI_Type_commit(&local_block);

    //----Rank 0 defines positions and counts of local blocks (2D-subdomains) on global matrix----//
    int * scatteroffset, * scattercounts;
    if (rank==0) {
        scatteroffset=(int*)malloc(size*sizeof(int));
        scattercounts=(int*)malloc(size*sizeof(int));
        for (i=0;i<grid[0];i++)
            for (j=0;j<grid[1];j++) {
                scattercounts[i*grid[1]+j]=1;
                scatteroffset[i*grid[1]+j]=(local[0]*local[1]*grid[1]*i+local[1]*j);
            }
    }

    // Debug print
    // if (rank==0) {
    //     printf("local[0]: %d, local[1]: %d\n", local[0], local[1]);
    //     for (i=0; i<global_padded[0]; i++) {
    //         for (j=0; j<global_padded[1]; j++)    
    //             printf("%f ", U[i][j]);
    //             //printf("%p ", (void *)&U[i][j]);
    //         printf("\n");
    //     }
    //     printf("\n");
    // }
    MPI_Barrier(CART_COMM);
    //----Rank 0 scatters the global matrix----//
    
    //----Rank 0 scatters the global matrix----//

	//*************TODO*******************//



	/*Fill your code here*/

	/*Make sure u_current and u_previous are
		both initialized*/
    zero2d(u_previous, local[0]+2, local[1]+2);
    zero2d(u_current, local[0]+2, local[1]+2);

    MPI_Barrier(CART_COMM);    
    MPI_Scatterv(&(U[0][0]), scattercounts, scatteroffset, global_block, &(u_current[1][1]), 1, local_block, 0, CART_COMM);
    
    //printf("Reached here 0\n");
    //************************************//

    // Copy u_current to u_previous to maintain initial border values
    for (i=0; i<local[0]+2; i++){
        for (j=0; j<local[1]+2; j++){
            u_previous[i][j]=u_current[i][j];
        }
    }

    if (rank==0)
        free2d(U);

     
	//----Define datatypes or allocate buffers for message passing----//
        
    //*************TODO*******************//

	/*Fill your code here*/

    double * right_send, * left_send, * up_send, * down_send, * right_receive, * left_receive, * up_receive, * down_receive;

	//************************************//

    //----Find the 4 neighbors with which a process exchanges messages----//

	//*************TODO*******************//
    int north, south, east, west;
    north = south = east = west = -1;
	/*Fill your code here*/
    
    // We are at rank_grid creds
    int north_g[2] = {rank_grid[0]-1, rank_grid[1]};
    if (rank_grid[0] > 0)
        MPI_Cart_rank(CART_COMM, north_g, &north);
    
    int south_g[2] = {rank_grid[0]+1, rank_grid[1]};
    if (rank_grid[0]+1 < grid[0])
        MPI_Cart_rank(CART_COMM, south_g, &south);

    int west_g[2] = {rank_grid[0], rank_grid[1]-1};
    if (rank_grid[1] > 0)
        MPI_Cart_rank(CART_COMM, west_g, &west);

    int east_g[2] = {rank_grid[0], rank_grid[1]+1};
    if (rank_grid[1]+1 < grid[1])
        MPI_Cart_rank(CART_COMM, east_g, &east);
    
	
    /*Make sure you handle non-existing
		neighbors appropriately*/
	//************************************//


    //---Define the iteration ranges per process-----//
	//*************TODO*******************//

	/*Fill your code here*/
    
    int i_min,i_max,j_min,j_max;
    j_min = ((west==-1)?1:0) + 1; // if we are at the first column -> i_min = 2
    i_min = ((north==-1)?1:0) + 1; // if we are at the first row -> j_min = 2

    if (south == -1) {
        i_max = ((global[0]%grid[0]!=0) ? (global[0]%grid[0]): local[0]);
    }
    else {
        i_max = local[0]+1;
    }
 
    if (east == -1) {
        j_max = (global[1]%grid[1]!=0) ? (global[1]%grid[1]): local[1];
    }
    else {
        j_max = local[1]+1;
    }
    
	/*Three types of ranges:
		-internal processes
		-boundary processes
		-boundary processes and padded global array
	*/

    up_send = (double *) calloc(local[1],sizeof(double));
    up_receive = (double *) calloc(local[1],sizeof(double));

    down_send = (double *) calloc(local[1],sizeof(double));
    down_receive = (double *) calloc(local[1],sizeof(double));

    right_send = (double *) calloc(local[0],sizeof(double));
    right_receive = (double *) calloc(local[0],sizeof(double));    
    
    left_send = (double *) calloc(local[0],sizeof(double));
    left_receive = (double *) calloc(local[0],sizeof(double));

	//************************************//

 	//----Computational core----//   
	gettimeofday(&tts, NULL);
    #ifdef TEST_CONV
    for (t=0;t<T && !global_converged;t++) {
    #endif
    #ifndef TEST_CONV
    #undef T
    #define T 256
    for (t=0;t<T;t++) {
    #endif
        //*************TODO*******************//

        /*Fill your code here*/

        /*Compute and Communicate*/

        /*Add appropriate timers for computation*/
        MPI_Request send_requests[4];
        MPI_Request receive_requests[4];



        
        // Copy data to send buffers
        for (i = 0; i < local[0]; i++) {
            right_send[i] = u_current[i+1][local[1]];
            left_send[i] = u_current[i+1][1];
        }        


        for (i = 0; i < local[1]; i++) {
            up_send[i] = u_current[1][i+1];
            down_send[i] = u_current[local[0]][i+1];
        }

        gettimeofday(&tes,NULL);
        int counter = 0;
        if (north != -1) {    
        //  Send north
            MPI_Isend(up_send,local[1],MPI_DOUBLE,north,north,CART_COMM, &send_requests[counter]);
            MPI_Irecv(up_receive,local[1],MPI_DOUBLE,north,rank,CART_COMM, &receive_requests[counter]);
            counter++;
        }

        if (south != -1) {
        //  Send south
            MPI_Isend(down_send,local[1],MPI_DOUBLE,south,south,CART_COMM, &send_requests[counter]);
            MPI_Irecv(down_receive,local[1],MPI_DOUBLE,south,rank,CART_COMM, &receive_requests[counter]);
            counter++;
        }

        if (west != -1) {
        //   Send west 
            MPI_Isend(left_send,local[0],MPI_DOUBLE,west,west,CART_COMM, &send_requests[counter]);
            MPI_Irecv(left_receive,local[0],MPI_DOUBLE,west,rank,CART_COMM, &receive_requests[counter]);
            counter++;
        }

        if (east != -1) {
        //   Send east 
            MPI_Isend(right_send,local[0],MPI_DOUBLE,east,east,CART_COMM, &send_requests[counter]);
            MPI_Irecv(right_receive,local[0],MPI_DOUBLE,east,rank,CART_COMM, &receive_requests[counter]);
            counter++;
        }

        MPI_Waitall(counter,receive_requests,MPI_STATUSES_IGNORE);               
        MPI_Waitall(counter,send_requests,MPI_STATUSES_IGNORE);               
        gettimeofday(&tef,NULL);
        tcomm+=(tef.tv_sec-tes.tv_sec)+(tef.tv_usec-tes.tv_usec)*0.000001;


        if (east!=-1)
            for (i = 0; i<local[0]; i++) {
                u_current[i+1][local[1]+1] = right_receive[i];
            }

        if (west!=-1)
            for (j = 0; j < local[0]; j++) {
                u_current[j+1][0] = left_receive[j];
            } 

        if (north!=-1)
            for (i = 0; i < local[1]; i++) {
                u_current[0][i+1]=up_receive[i];
            }
        
        if (south!=-1)
            for (j = 0; j < local[1]; j++) {
                u_current[local[0]+1][j+1] = down_receive[j];
            }

        // MPI_Barrier(CART_COMM);

		swap=u_previous;
		u_previous=u_current;
		u_current=swap;

        // Computation timer
		gettimeofday(&tcs,NULL);
        for (i=i_min; i<i_max; i++)
            for (j=j_min; j<j_max; j++) {
                u_current[i][j]=(u_previous[i-1][j]+u_previous[i+1][j]+u_previous[i][j-1]+u_previous[i][j+1])/4.0;
            }
        
        gettimeofday(&tcf,NULL);
        tcomp+=(tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;
        

        //************************************//


        // MPI_Barrier(CART_COMM);
        #ifdef TEST_CONV
        if (t%C==0) {
            gettimeofday(&tss,NULL);
            //*************TODO**************//
            /*Test convergence*/
            // Local conv check, converged = 1 if conv else converged = 0
            converged = converge(u_previous, u_current, i_min, i_max-1, j_min, j_max-1);
            MPI_Allreduce(&converged, &global_converged, 1, MPI_INT, MPI_MIN, CART_COMM);
            gettimeofday(&tsf,NULL);
            tconv+=(tsf.tv_sec-tss.tv_sec)+(tsf.tv_usec-tss.tv_usec)*0.000001;

        }
        #endif
        // Just for debug on no conv test 
        // Error still exists even with barrier here
        // MPI_Barrier(CART_COMM);
    }

    free(up_send);
    free(up_receive);

    free(down_send);
    free(down_receive);

    free(right_send);
    free(right_receive);    
    
    free(left_send);
    free(left_receive);

    MPI_Barrier(CART_COMM);
    gettimeofday(&ttf,NULL);

    ttotal=(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001;

    MPI_Reduce(&ttotal,&total_time,1,MPI_DOUBLE,MPI_MAX,0,CART_COMM);
    MPI_Reduce(&tcomp,&comp_time,1,MPI_DOUBLE,MPI_MAX,0,CART_COMM);
    MPI_Reduce(&tcomm,&comm_time,1,MPI_DOUBLE,MPI_MAX,0,CART_COMM);
    // #ifdef TEST_CONV
    //     MPI_Reduce(&tconv,&conv_time,1,MPI_DOUBLE,MPI_MAX,0,CART_COMM);
    // #endif


    //----Rank 0 gathers local matrices back to the global matrix----//

    if (rank==0) {
        U=allocate2d(global_padded[0],global_padded[1]);
    }


    //*************TODO*******************//
    /*Fill your code here*/
    // Gather
    MPI_Gatherv(&u_current[1][1], 1, local_block, &U[0][0], scattercounts, scatteroffset, global_block, 0, CART_COMM);
    
    //************************************//


    //----Printing results----//

    //**************TODO: Change "Jacobi" to "GaussSeidelSOR" or "RedBlackSOR" for appropriate printing****************//
    if (rank==0) {
        // #ifdef TEST_CONV
        //     printf("Jacobi X %d Y %d Px %d Py %d Iter %d ComputationTime %lf TotalTime %lf CommunicationTime %lf ConvTime %lf midpoint %lf\n",global[0],global[1],grid[0],grid[1],t,comp_time,total_time,comm_time,conv_time,U[global[0]/2][global[1]/2]);
        // #else
            printf("Jacobi X %d Y %d Px %d Py %d Iter %d ComputationTime %lf TotalTime %lf CommunicationTime %lf midpoint %lf\n",global[0],global[1],grid[0],grid[1],t,comp_time,total_time,comm_time,U[global[0]/2][global[1]/2]);
        // #endif
        #ifdef PRINT_RESULTS
            char * s=malloc(50*sizeof(char));
            sprintf(s,"TEST_CONV_resJacobiMPI_%dx%d_%dx%d_iters_%d",global[0],global[1],grid[0],grid[1],t);
            fprint2d(s,U,global[0],global[1]);
            free(s);
        #endif

    }
    
    MPI_Finalize();
    return 0;    
}
