// Parallelise your program

/* each process has unique ID --> rank
 each works on different part of the problem


 link with MPI libraries or compile with mpicxx

 run with mpiexec command

	mpicxx hello.cpp -o hello
	mpiexcec -np 2 ./hello

	

each message contains:
- sender rank
- destination rank
- message datatype
- message size
- message location
- message tag

few libraries needed for scalapack:
mpicxx myprog
.
cpp
-
lscalapack
-
openmpi
-
lblacs
-
openmpi
\
-
lblacsCinit
-
openmpi
-
llapack
-
lblas


 */



void T4_inputs(){
	std::cout << " ------------ TASK 4 -------------" << std::endl;
	int size;
    int rank;                   
    size = MPI_Comm_size(MPI_COMM_WORLD, &rank);
    std::cout <<"Number of processes: " << size << std::endl;		
    rank = MPI_Comm_rank(MPI_COMM_WORLD, &size);


	return;
}


