## Semester Exercises in Parallel Processing Systems  
*National Technical University of Athens* 

> **Disclaimer**: A significant portion of the code was provided, and our primary task was to make additions or modifications to specific sections, guided by `//TODO` comments throughout the codebase.

#### Team 38
| AM      | Name                    |
|---------|-------------------------|
| el20034 | Μαρντιροσιάν Φίλιππος    |
| el20874 | Μπέλσης Παναγιώτης      |
| el20109 | Θεοδοσίου Γεώργιος      |

### Exercise 1: Conway's Game of Life  
Parallelized the cellular automaton using OpenMP on shared-memory systems.

### Exercise 2: K-means & Floyd-Warshall  
Optimized K-means clustering with OpenMP (shared/copied clusters) and pthread locks (mutexes, spinlocks). Parallelized the recursive Floyd-Warshall algorithm using OpenMP tasks, testing scalability for different matrix sizes.

### Exercise 3: GPU K-means with CUDA  
Implemented and optimized K-means on NVIDIA GPUs (naive, transpose, shared memory versions). Analyzed multiple block sizes and offload efficiency, comparing GPU execution, data transfers, and CPU compute times.

### Exercise 4: MPI K-means & Heat Diffusion  
Parallelized K-means and heat diffusion solvers (Jacobi, Gauss-Seidel SOR) using MPI. Tested scalability for different numbers of processes on distributed grids, focusing on computation vs. communication trade-offs.  

## Technologies Used  
OpenMP, CUDA, MPI, Pthreads, C/C++, Bash Scripts (.sh)  
