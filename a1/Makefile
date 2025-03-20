all: omp_gol

omp_gol: Game_Of_Life.c
	gcc -O3 -fopenmp -o GameOfLife Game_Of_Life.c

clean: 
	rm GameOfLife
	rm performance_results.txt
