CC=mpicc
CFLAGS=-O3 -lm -Wall
RES=-DPRINT_RESULTS
CONV=-DTEST_CONV

all: jacobi

jacobi: mpi_jacobi.c utils.c
	# $(CC) $(CFLAGS) $(RES) $(CONV) mpi_jacobi.c utils.c -o jacobi
	# $(CC) $(CFLAGS) $(CONV) mpi_jacobi.c utils.c -o jacobi
	$(CC) $(CFLAGS) mpi_jacobi.c utils.c -o jacobi
clean:
	rm jacobi
