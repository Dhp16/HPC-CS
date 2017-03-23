/*  
 Imperial College London
 HPC Assignment Task 5
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 4

 The main loop is very close to that designed for task 2, with the addition of communication between the processes

 Inspiration from example pbsv.cpp provided
*/


/* Function: Build_PARA_K builds a stiffness matrix of the correct size for each process and be passed to pdgbsv*/
void Build_PARA_K(double* PARA_K, const int PARA_nodes, const int Nx, const double A, const double E, const double I, const double l){	
	double Ke[36] = {0};
	pp_Ke(Ke, A, E, I, l);
	const int k = 8; 								// rows of 0s above
	const int n = 3 ;
	double* PARA_Ke_W0s = new double[(PARA_nodes+4)*n*(9+k)];		// with the banded 0s
	Global_Stiffness_Matrix(PARA_Ke_W0s, Ke, PARA_nodes+4, k);
	Matrix_Transformer_Imp(PARA_Ke_W0s, PARA_nodes+4, k);

	int counter = 0;
	int i = 102;
	while(counter < PARA_nodes*n*(9+k)){
		PARA_K[counter] = PARA_Ke_W0s[i];
		counter++;
		i++;
	}

	//delete[] PARA_Ke_W0s;

	return;
}


void T5_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){
	std::cout << " ------------ TASK 5 ------------- " << std::endl;
 	
 	const int n = 3;											// columns per node
	const double beta = 0.250; 									// constant beta
	const double gamma = 0.5;									// constant
	double delta_t = T/Nt;										// time step length (s)

	double Fe[6] = {0};											// elemental force vector
	pp_Fe(Fe, l, Nx, Qy, Qx);									// populate the elemental force vector

	int rank;
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size); 						// access the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);						// access the rank of a process

 	//const int PS   = Nx;			// Problem size
 	//const int nb = PS/size;

 	const int midinc_rank = size/2;

  	const int PS = (Nx+1)*3;			// Problem size
 	const int nb = (Nx+1)*3/size;

	//const int nb   = PS/size;		// Block PARA_size --> need to make sure its a multiple of 3
    const int kl   = 4;   			// Number of lower diagonals
    const int ku   = 4;   			// Number of upper diagonals
    const int lda  = 1 + 2*kl + 2*ku; // Leading dimension (num of rows)
    const int ldb  = nb;
    const int nrhs = 1;
    const int ja   = 1;             // Offset index in global array (col)
    const int ib   = 1;             // Offset index in global array (row)
    const int lwork= (nb+ku)*(kl+ku)+6*(kl+ku)*(kl+2*ku) + std::max(nrhs*(nb+2*kl+4*ku), 1);
    int info = 0;	

    int nrow = 4;		// 1 (in example) 
    int ncol = 4;		// 2 (in example)
    char order = 'R';
	int ctx;
    int mype;
    int npe;
    int myrow;
    int mycol;
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get( 0, 0, &ctx );
    Cblacs_gridinit( &ctx, &order, 1, npe );
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol); 

     // Create descriptors for matrix and RHS vector storage			// yet to be sorted out
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
    descb[4] = 0;           // Process row/column						// ----------> what!?!?!?
    descb[5] = nb;          // Local size
    descb[6] = 0;           // Reserved

    int PARA_nodes = nb/3;								 // includes the 0s needed 								// -----> NEED TO RESET
    int PARA_size = PARA_nodes*3;						 // A equivalent --> need to build it
    double* PARA_SolBlock = new double[nb*3];

    if(rank == 0){
    	std::cout << "Nx: " << Nx << std::endl;
    	std::cout << "PS: " << PS << std::endl;
    	std::cout << "PARA_node: " << PARA_nodes << std::endl;
    	std::cout << "nb: " << nb << std::endl;
    	std::cout << "midinc_rank" << midinc_rank << std::endl;
    }
    // Specify the force with concentrated for only one force

    // -------------------------------------- Equation 13 ----------------------------------------
    double* K_effective = new double[PARA_nodes*3*lda];				 // A equivalent --> need to build it     /// wrong size 
    Build_PARA_K(K_effective, nb, Nx, A, E, I, l);
    double* G_Mm = new double[PARA_size];
	Global_Mass_Matrix(G_Mm, PARA_nodes, rho, A, l);
	K_Effective(K_effective, G_Mm, PARA_nodes, delta_t);

    // -------------------------------------- Setting up the LOOP --------------------------------
	double* DIS = new double[PARA_size];
	double* VEL = new double[PARA_size];
	double* ACC = new double[PARA_size];
	initialise_dynamic_array(DIS, PARA_nodes, n, 1);
	initialise_dynamic_array(VEL, PARA_nodes, n, 1);
	initialise_dynamic_array(ACC, PARA_nodes, n, 1);
	double* DIS_plus = new double[PARA_size];
	double* VEL_plus = new double[PARA_size];
	double* ACC_plus = new double[PARA_size];
	initialise_dynamic_array(DIS_plus, PARA_nodes, n, 1);
	initialise_dynamic_array(VEL_plus, PARA_nodes, n, 1);
	initialise_dynamic_array(ACC_plus, PARA_nodes, n, 1);

    // -------------------------------------- Running the LOOP  --------------------------------
	for(int i = 0; i < Nt; i++){
		double t_now = i*delta_t;
		double* Fi = new double[PARA_size];
		initialise_dynamic_array(Fi, PARA_nodes, n, 1);
		if(midinc_rank != 0){
			Global_Force_Vector(Fi, Fe, PARA_nodes, 0, 0, (t_now+delta_t)/T);
			if(rank == midinc_rank)
				Fi[PARA_size-2] = Fi[PARA_size-2] + Fy*(t_now+delta_t)/T;							// add the concentrated load to the central node 
		}
		else{
			Global_Force_Vector(Fi, Fe, PARA_nodes, 0, Fy, (t_now+delta_t)/T);
		}

		// ---------------------- Equation 12 -------------------
		double* A = new double[PARA_size];
		for(int i = 0; i < PARA_size; i++)
			A[i] = (1/(beta*pow(delta_t, 2))) * DIS[i] + (1/(beta*delta_t)) * VEL[i] + (0.5/beta-1)*ACC[i];

		double* B = new double[PARA_size];
		double* B_copy = new double[PARA_size];
		for(int i = 0; i < PARA_size; i++)
			B[i] = Fi[i] + A[i]*G_Mm[i];

		// F77NAME(dgbsv) (vec_size, 4, 4, 1, K_effective_copy, 1+2*4+4, ipiv, B, vec_size, &info);	// result should be ported to DIS_plus10--> save KEff for next loop
 
		double* K_effective_copy = new double[PARA_size*lda];
		Matrix_Copy(K_effective_copy, K_effective, PARA_size*lda);
		double* wk   = new double[lwork];
		int* ipiv = new int[nb];
    	// Solve the system A * y' = y (i.e. RHS vector replaced by solution)
    	Matrix_Copy(B_copy, B, PARA_size);
    	F77NAME(pdgbsv) (nb, kl, ku, nrhs, K_effective, ja, desca, ipiv, B, ib, descb, wk, lwork, &info);
    	delete[] K_effective_copy;

		// ---------------------- Equation 10 ------------------- 
    	for(int i = 0; i < PARA_size; i++)
    		ACC_plus[i] = ((1/(beta*pow(delta_t, 2)))*(DIS_plus[i] - DIS[i])) - (1/(beta*delta_t))*VEL[i] - ((1/(2*beta)-1)*ACC[i]);   

    	// ----------------------- Equation 11 -------------------
    	for(int i = 0; i < PARA_size; i++)
    		VEL_plus[i] =       VEL[i]      +         (delta_t*(1-gamma))*ACC[i]     +    delta_t*gamma*ACC_plus[i];

    	// ---------------------- Assigning for next loop --------------
    	Matrix_Copy(DIS, DIS_plus, PARA_size);
    	Matrix_Copy(VEL, VEL_plus, PARA_size);
    	Matrix_Copy(ACC, ACC_plus, PARA_size);
    	initialise_dynamic_array(DIS_plus, PARA_nodes, n, 1);	// making sure these are returned to 0 before next iteration
		initialise_dynamic_array(VEL_plus, PARA_nodes, n, 1);
		initialise_dynamic_array(ACC_plus, PARA_nodes, n, 1);

		if(rank == 0){
			std::cout << "\n--------------------- Time step: " << i << " ---------------" << std::endl;
			for(int i = 1; i < PARA_size; i+= 3)
				std::cout << std::setw(10) << Fi[i];
			std::cout << "\n" << info << " | ";
			for(int i = 1; i < PARA_size; i+= 3)
				std::cout << std::setw(10) << DIS[i];
			std::cout << "\n-------------------------------------" << std::endl;
		}

		delete[] A;
		delete[] B;
		delete[] B_copy;
		delete[] ipiv;
		delete[] wk;
    }


	if(rank == 0) std::cout << "END OF LOOP" << std::endl;

  	// use send receive to get the displacement from each process and collate into ONE!

    Cblacs_gridexit( ctx );


    if(rank == 0) std::cout << "DIS" << std::endl;delete[] DIS;
    if(rank == 0) std::cout << "VEL" << std::endl;delete[] VEL;
    if(rank == 0) std::cout << "ACC" << std::endl;delete[] ACC;
    if(rank == 0) std::cout << "DIS_plus" << std::endl;delete[] DIS_plus;
    if(rank == 0) std::cout << "VEL_plus" << std::endl;delete[] VEL_plus;
    if(rank == 0) std::cout << "ACC_plus" << std::endl;delete[] ACC_plus;

    if(rank == 0) std::cout << "K_effective" << std::endl;delete[] K_effective;				 // A equivalent --> need to build it
    if(rank == 0) std::cout << "G_Mm" << std::endl;delete[] G_Mm;
    if(rank == 0) std::cout << "PARA_SolBlock" << std::endl;delete[] PARA_SolBlock;

    return;
}

// force always in the next to last node in P/2 - 1


