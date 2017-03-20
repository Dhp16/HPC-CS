/* 
 Imperial College London
 HPC Assignment Task 3
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 4

 The main loop is very close to that designed for task 2, with the addition of communication between the processes
*/



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

 */

//void PARA_build_M(const double rho, const double A, const double l, const int PARA_size, const int rank){

// /} 
void PARA_Force_Vector(double* PARA_G_Fe, const double l, const int PARA_Nx, const int Nx, const double Fy, const double Qy, const double Qx, const int rank,
	const double scalar){
	double* G_Fe = new double[Nx*3];
	double Fe[6] = {0};
    pp_Fe(Fe, l, Nx, Qy, Qx);
   	Global_Force_Vector(G_Fe, Fe, Nx, 0, Fy, scalar);
   	if(rank == 0){
   		for(int i = 0; i < PARA_Nx*3; i++)
   			PARA_G_Fe[i] = G_Fe[i];
   	}
   	else if(rank == 1){
   		int counter = 0;
   		for(int i = PARA_Nx*3-5*3; i < Nx*3; i++){
    		PARA_G_Fe[counter] = G_Fe[i];
    		counter++;
    	}
    }
    delete[] G_Fe;
}    	


void T4_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){
	std::cout << " ------------ TASK 4 -------------" << std::endl;
	int size;
    int rank;
    int PARA_Nx = Nx/2+3;  
    int PARA_size = 3*PARA_Nx; 
    int n = 3;
    double delta_t = T/Nt;

    std::cout << "PARA_size: " << PARA_size << std::endl;                 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout <<"Number of processes: " << size << std::endl;		
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	std::cout << "Processor: " << rank << std::endl;

    double* PARA_G_Mm = new double[PARA_size];
    Global_Mass_Matrix(PARA_G_Mm, PARA_Nx, rho, A, l);
    
	double* PARA_G_Ke = new double[PARA_size*9];
	double Ke[36] = {0};				
	pp_Ke(Ke, A, E, I, l);	
	Global_Stiffness_Matrix(PARA_G_Ke, Ke, PARA_Nx, 0 );

    double Fe[6] = {0};
    pp_Fe(Fe, l, PARA_Nx, Qy, Qx);


    // Initialising constants in the integration scheme	-------------------------------------------------------------------------------------------
    double* COEF1 = new double[PARA_size];			// COEF1 = delta_T^2/[M]
	double* COEF2 = new double[PARA_size*9];			// COEF2 = [K] - 2/dealtaT^2 [M]
	double* COEF3 = new double[PARA_size];			// COEF3 = -1/deltat^2 * [M] 

	Build_Coef1(COEF1, PARA_G_Mm, PARA_Nx, delta_t);
	Build_Coef2(COEF2, PARA_G_Ke, PARA_G_Mm, PARA_Nx, delta_t);        // at this stage COEF2 = [K] - 2[M]/delta_t^2  -    banded - column Major
	Build_Coef3(COEF3, PARA_G_Mm, PARA_Nx, delta_t);
	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Initialising for loop ----------------------------------------------------------------------------------------------------------------------
	double* U_minus = new double[PARA_size];
	double* U_now = new double[PARA_size];
	double* U_plus = new double[PARA_size];  
	initialise_dynamic_array(U_minus, PARA_Nx, n, 1);
	initialise_dynamic_array(U_now, PARA_Nx, n, 1);
	initialise_dynamic_array(U_plus, PARA_Nx, n, 1); 
	// --------------------------------------------------------------------------------------------------------------------------------------------

	std::cout << "Process " << rank << " starting loop." << std::endl;

 	for(int i = 0; i < Nt; i++){
 		double t_now = i*delta_t;
 		double* PARA_G_Fe = new double[PARA_size];
 		PARA_Force_Vector(PARA_G_Fe, l, PARA_Nx, Nx, Fy, Qy, Qx, rank, t_now/T);

		double* Y2 = new double[PARA_size];
		double* Y3 = new double[PARA_size]; 
		initialise_dynamic_array(Y2, PARA_Nx, n, 1);
		initialise_dynamic_array(Y3, PARA_Nx, n, 1);

		F77NAME(dgbmv)('n', PARA_size, PARA_size, 4, 4, 1.0, COEF2, 9, U_now, 1, 0.0, Y2, 1);		// Y2 = COEF2*U_n

		Diagonal_by_Vector(Y3, COEF3, U_minus, PARA_size);	
		double* result = new double[PARA_size];
		initialise_dynamic_array(result, PARA_Nx, n, 1);
		adding_three_arrays(result, PARA_G_Fe, Y2, Y3, PARA_size);

		Diagonal_by_Vector(U_plus, COEF1, result, PARA_size);

		Matrix_Copy(U_minus, U_now, PARA_size);
		Matrix_Copy(U_now, U_plus, PARA_size);
  
		// collate data and exchange ---------------------------------------
     MPI_Barrier(MPI_COMM_WORLD);
   		if( i < Nt){
   			//if(rank == 1){
   			//	std::cout << "\n-------------- Time step: " << i << " Tnow: " << t_now  << "-----------";	
        //}
   			
        double Shared[12];

   			if(rank == 0){
            int counter = 0;
    		    for(int i = PARA_size-12; i < PARA_size; i++){
    		    	Shared[counter] = U_now[i];
              counter++;
            } 
    		}
    		else if (rank == 1){
    		    for(int i = 0; i < 12; i++)
    		    	Shared[i] = U_now[i];
    		}
	
    		double post_box[12];
    		if(rank == 0)
    		    MPI_Send(&Shared, 12, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    		else if(rank == 1){
    		    MPI_Send(&Shared, 12, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    		}
    
        MPI_Barrier(MPI_COMM_WORLD);		
         
    		if(rank == 0){
    		    MPI_Recv(&post_box, 12, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            int counter = 0;
            for(int i = PARA_size-6; i < PARA_size; i++){
              U_now[i] = post_box[counter];
              counter++;
            }
    		}
    		else if(rank == 1){
    		    MPI_Recv(&post_box, 12, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int i = 0; i < 6; i++)
              U_now[i] = post_box[i];
    		} 
    	}
		// -----------------------------------------------------------------

		delete[] result;
		delete[] Y2;
		delete[] Y3;
 		delete[] PARA_G_Fe;
 	}
 
 	if(rank == 1){
 		std::cout << "\n\n\n\n\n\nRAnk: " << rank << " solution" << std::endl;
 		for(int i = 1; i < PARA_size; i+=3)
 			std::cout << U_now[i] << "   ";
 	}

 	// use process 0 to write results to files

	//delete[] PARA_G_Mm
 	delete[] U_minus;
 	delete[] U_now;
 	delete[] U_plus;

    delete[] PARA_G_Mm;
    delete[] PARA_G_Ke; 
	return;
}


// at the very end 

// wait for next processs at the end of each loop to pass the information