/*  
 Imperial College London
 HPC Assignment Task 5
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 4

 The main loop is very close to that designed for task 2, with the addition of communication between the processes

 Inspiration from example pbsv.cpp provided
*/



# include "task5.h"
# include "tools.h"
# include "task1.h"
# include "task2.h"
# include <vector>
# include <string>
# include <mpi.h>  

 
void T5_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){
	
	std::cout << "\n---------------------- Task 5 ------------------------" << std::endl;

	/*
		Once the boundaries are taken off there is an uneven number of nodes to distribute evenly across a pair number of processes.
		In order to maintain an equal block size across processes, the last process is allocated one more node that needed. 
		This node is populated with 1s in its diagonal and 0s elsewhere so as not to be ill-defined.
	*/

	const double gamma = 0.5;					// constants
	const double beta = 0.25;					
	const int n = 3;
	const int k = 8;							// number of lines of 0s needed by pdgbsv
	const double delta_t = T/Nt;				// delta_t

	// Populate 
	double Fe[6] = {0};
	pp_Fe(Fe, l, Nx, Qy, Qx);
	double Ke[36] = {0};
	pp_Ke(Ke, A, E, I, l);

	int rank;									
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size); 		// access the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int PARA_Nx;
    int midinc_rank;					// identify the rank that will include the concentrated load

    if(size > 1){						// set the size of the arrays needed for 
    	PARA_Nx = (Nx+1)/size;
		midinc_rank = size/2-1;
	}
	else
		PARA_Nx = Nx-1;

    const int PS   = Nx*3;				// Problem size
	const int nb   = Nx*3/size;			// Block size
    const int kl   = 4;   			    // Number of lower diagonals
    const int ku   = 4;   			    // Number of upper diagonals
    const int lda  = 1 + 2*kl + 2*ku;   // Leading dimension (num of rows)
    const int ldb  = nb;
    const int nrhs = 1;
    const int ja   = 1;            		// Offset index in global array (col)
    const int ib   = 1;            		// Offset index in global array (row)
    const int lwork= (nb+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku) + std::max(nrhs*(nb+2*kl+4*ku), 1);
    int info = 0;

    // Set up the CBLACS context
    int nrow = 1;
    int ncol = size;
    char order = 'C';
	int ctx;
    int mype;
    int npe;
    int myrow;
    int mycol;
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get( 0, 0, &ctx );
    Cblacs_gridinit( &ctx, &order, 1, npe );
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);

    // Create descriptors for matrix and RHS vector storage
    int desca[7]; 
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = ctx;         // Context
    desca[2] = PS;          // Problem size
    desca[3] = nb;          // Blocking
    desca[4] = 0;           // Process row/column
    desca[5] = lda;         // Local size
    desca[6] = 0;           // Reserved

    int descb[7]; 
    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = ctx;         // Context
    descb[2] = PS;          // Problem size
    descb[3] = nb;          // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = nb;          // Local size
    descb[6] = 0;           // Reserved
    
    // -------------------------------------- Equation 13 ----------------------------------------
	int counter = 0;
    // K matrix for each process ------------
	double* PARA_Ke_0s = new double[(PARA_Nx+4)*n*(9+k)];
	Global_Stiffness_Matrix(PARA_Ke_0s, Ke, PARA_Nx+4, k);
	Matrix_Transformer_Imp(PARA_Ke_0s, PARA_Nx+4, k);
	double* PARA_Ke = new double[nb*(9+k)];
	int i = 102;
	while(counter < nb*(9+k)){
		PARA_Ke[counter] = PARA_Ke_0s[i];
		counter++;
		i++;
	}
	delete[] PARA_Ke_0s;

	// Adding the Mass factor ---------------
	double* G_Mm = new double[nb];
	Global_Mass_Matrix(G_Mm, PARA_Nx, rho, A, l);
	if(size != 1 && rank == size-1){					
		G_Mm[nb-3] = 0;								// setting 0s in excess nodes 
		G_Mm[nb-2] = 0;
		G_Mm[nb-1] = 0;
	}

	counter = 0;
	for(int i = 0; i < nb*(9+k); i++){
		if((i+5) % 17 == 0){
			PARA_Ke[i] = PARA_Ke[i] + 1/(beta*pow(delta_t, 2))*G_Mm[counter];
			counter++;
		}
	}

	// -------------------------------------	Setting the 1s and 0s to the diagonal on excess Ke -----------
	if(size != 1 && rank == size-1){
		for(int i = (PARA_Nx-1)*3*17; i < nb*17; i++){
			if((i+5) % 17 == 0){
				PARA_Ke[i] = 1.0;
			}
			else
				PARA_Ke[i] = 0.0;
		}
	}
    // -------------------------------------- Setting up the LOOP --------------------------------
	double* DIS = new double[nb];									// initialse 
	double* VEL = new double[nb];
	double* ACC = new double[nb];
	initialise_dynamic_array(DIS, nb);
	initialise_dynamic_array(VEL, nb);
	initialise_dynamic_array(ACC, nb);
	double* DIS_plus = new double[nb];
	double* VEL_plus = new double[nb];
	double* ACC_plus = new double[nb];
	initialise_dynamic_array(DIS_plus, nb);
	initialise_dynamic_array(VEL_plus, nb);
	initialise_dynamic_array(ACC_plus, nb);

	// -------------------------------------- Time loop ----------------------------------
	for(int i = 0; i < Nt; i++){
		double t_now = i*delta_t;

		// ---- Force --------------------------------------------------------------------- 
		double* Fi = new double[nb];
		initialise_dynamic_array(Fi, nb);
		if(size > 1){
			Global_Force_Vector(Fi, Fe, PARA_Nx, 0, 0, (t_now+delta_t)/T);	// do not include the concetrated force when populating the force for each process
			if(rank == midinc_rank)
				Fi[nb-2] += Fy*(t_now+delta_t)/T;							// add the concentrated load to the central node 
			else if(rank == size -1){
				Fi[nb-3] = 0;
				Fi[nb-2] = 0;
				Fi[nb-1] = 0;
			}
		}
		else if(size == 1){
			Global_Force_Vector(Fi, Fe, PARA_Nx, 0, Fy, (t_now+delta_t)/T);	// Force for the single process is the same as the global force
		}

		// ----------------------- RHS of Equation 12 ----------------------------------
		double* A = new double[nb];
		for(int i = 0; i < nb; i++)
			A[i] = (1/(beta*pow(delta_t, 2))) * DIS[i] + (1/(beta*delta_t)) * VEL[i] + (0.5/beta-1)*ACC[i];

		double* B = new double[nb];
		for(int i = 0; i < nb; i++)
			B[i] = Fi[i] + A[i]*G_Mm[i];			// adding the force

		// clear after equation
		delete[] A;

		// Setting up for solvernb
		double* PARA_Ke_copy = new double[nb*(9+k)];						// the solver modifies the 4th entry, hence it is copied before being passed
		Matrix_Copy(PARA_Ke_copy, PARA_Ke, nb*(9+k));
		int* ipiv = new int[PS];
		double* wk = new double[lwork]; 

		// ....................... Using Solver .................
		F77NAME(pdgbsv) (PS, kl, ku, nrhs, PARA_Ke_copy, ja, desca, ipiv, B, ib, descb, wk, lwork, &info);

		// Assign B to solution
		Matrix_Copy(DIS_plus, B, nb);

		// Clear after solver
		delete[] B;
		delete[] PARA_Ke_copy;
		delete[] ipiv;	// 
		delete[] wk;	//

		// ---------------------- Equation 10 ------------------- 
    	for(int i = 0; i < nb; i++)
    		ACC_plus[i] = ((1/(beta*pow(delta_t, 2)))*(DIS_plus[i] - DIS[i])) - (1/(beta*delta_t))*VEL[i] - ((1/(2*beta)-1)*ACC[i]);   

    	// ----------------------- Equation 11 -------------------
    	for(int i = 0; i < nb; i++)
    		VEL_plus[i] = VEL[i]      +         delta_t*(1-gamma)*ACC[i]     +    delta_t*gamma*ACC_plus[i];

    	Matrix_Copy(DIS, DIS_plus, nb);
    	Matrix_Copy(VEL, VEL_plus, nb);
    	Matrix_Copy(ACC, ACC_plus, nb);
    	initialise_dynamic_array(DIS_plus, nb);	// making sure these are returned to 0 before next iteration
		initialise_dynamic_array(VEL_plus, nb);
		initialise_dynamic_array(ACC_plus, nb);

		// clear memory
		delete[] Fi;

		// clear for Output
	} 

	MPI_Barrier(MPI_COMM_WORLD);
	// Clearing memory
	delete[] G_Mm;
	delete[] VEL;
	delete[] ACC;
	delete[] DIS_plus;
	delete[] VEL_plus;
	delete[] ACC_plus;
	delete[] PARA_Ke;


	/*
		If the task is running on one processor only, the latter is then passed on to the following
	*/

	if(size == 1)
		write_task_solution(DIS, nb, "5_size1");
	else{
		if(size == 2){																						// Solution sorting specific to 2 processes
			if(rank == 1){
				double* Parcel = new double[nb-3];
				Matrix_Copy(Parcel, DIS, nb-3);
				MPI_Send(Parcel, nb-3, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
			}
			if(rank == 0){
				double* Post_Box = new double[nb-3];
				MPI_Recv(Post_Box, nb-3, MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				
				double* Solution = new double[nb*size-3];
				for(int i = 0; i < nb; i++)
					Solution[i] = DIS[i];

				int counter = 0;
				for(int i = nb; i < nb*size-3; i++){
					Solution[i] = Post_Box[counter];
					counter++;
				}
				write_task_solution(Solution, nb*2-3, "5_size2");				
			}				
		}
		else if(size > 2){																					// Solution sorting for 4 or mor processes
			if(rank != 0 ){
				double* Parcel = new double[nb];
				Matrix_Copy(Parcel, DIS, nb);
				MPI_Send(Parcel, nb, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			if(rank == 0){
				double* Post_Box = new double[nb];
				double* Solution = new double[nb*size];
				int j = 0;
				for(; j < nb; j++)
					Solution[j] = DIS[j];

				for(int i = 1; i < size; i++){
					MPI_Recv(Post_Box, nb, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(int i = 0; i < nb; i++){
						Solution[j] = Post_Box[i];
						j++;
					}
				}
				write_task_solution(Solution, nb*size-3, "5_size4");				
			}
						
		}
	}

	// Clearing memory
	delete[] DIS;

    Cblacs_gridexit( ctx );

    if(rank == 0)
    	std::cout << "\n------------------- End of Task 5 ----------------------" << std::endl;

    return;
}

// force always in the next to last node in P/2 - 1


