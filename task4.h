/* 
 Imperial College London
 HPC Assignment Task 4
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

void write_task4(const int Nx, double* Solution, const int vec_size, const double l){
  std::ofstream file;
  std::string file_name = "Task4_data.txt"; 
  file.open(file_name);

  file << 0 << std::endl;
  for(int i = 1; i < vec_size; i+=3)
    file << Solution[i] << std::endl; 
  file << 0 << std::endl;

  for(int i = 0; i < Nx+2; i++)
    file << l*i << std::endl;

  file.close();

  std::cout << "\n ---------------- Written Task 4 results to file: " << file_name <<". ---------------------"<< std::endl;
  return;
}

void PARA_Force_Vector(double* PARA_G_Fe, const double l, const int PARA_Nx, const int Nx, const double Fy, const double Qy, const double Qx, const int rank, const double scalar){
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
  std::cout << "---------------------- Task 4 ------------------------" << std::endl;
	const bool output = false;
  int size;
  int rank;
  
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int vec_size = Nx*3;
  int PARA_Nx;
  if(size > 1)
    PARA_Nx = (Nx+1)/2+2;  
  else
    PARA_Nx = Nx;
  int PARA_size = 3*PARA_Nx; 
  int n = 3;
  double delta_t = T/Nt;                 

  if(output){
    std::cout <<"PARA_size: " << PARA_size << std::endl;   
    std::cout <<"Number of processes: " << size << std::endl;
  }

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
	double* COEF2 = new double[PARA_size*9];		// COEF2 = [K] - 2/dealtaT^2 [M]
	double* COEF3 = new double[PARA_size];			// COEF3 = -1/deltat^2 * [M] 

	Build_Coef1(COEF1, PARA_G_Mm, PARA_Nx, delta_t);
	Build_Coef2(COEF2, PARA_G_Ke, PARA_G_Mm, PARA_Nx, delta_t);        // at this stage COEF2 = [K] - 2[M]/delta_t^2  -    banded - column Major
	Build_Coef3(COEF3, PARA_G_Mm, PARA_Nx, delta_t);
	// --------------------------------------------------------------------------------------------------------------------------------------------

	// Initialising for loop ----------------------------------------------------------------------------------------------------------------------
	double* U_minus = new double[PARA_size];
	double* U_now = new double[PARA_size];
	double* U_plus = new double[PARA_size];  
	initialise_dynamic_array(U_minus, PARA_size);
	initialise_dynamic_array(U_now, PARA_size);
	initialise_dynamic_array(U_plus, PARA_size); 
	// --------------------------------------------------------------------------------------------------------------------------------------------

  if(output)
	 std::cout << "Process " << rank << " starting loop." << std::endl;

 	for(int i = 0; i < Nt; i++){
 		double t_now = i*delta_t;
 		double* PARA_G_Fe = new double[PARA_size];
 		if(size != 1)
      PARA_Force_Vector(PARA_G_Fe, l, PARA_Nx, Nx, Fy, Qy, Qx, rank, t_now/T);
    else
      Global_Force_Vector(PARA_G_Fe, Fe, PARA_Nx, 0, Fy, t_now/T);

		double* Y2 = new double[PARA_size]; 
		double* Y3 = new double[PARA_size]; 
		initialise_dynamic_array(Y2, PARA_size);
		initialise_dynamic_array(Y3, PARA_size);

		F77NAME(dgbmv)('n', PARA_size, PARA_size, 4, 4, 1.0, COEF2, 9, U_now, 1, 0.0, Y2, 1);		// Y2 = COEF2*U_n

		Diagonal_by_Vector(Y3, COEF3, U_minus, PARA_size);	
		double* result = new double[PARA_size];
		initialise_dynamic_array(result, PARA_size);
		adding_three_arrays(result, PARA_G_Fe, Y2, Y3, PARA_size);

		Diagonal_by_Vector(U_plus, COEF1, result, PARA_size);
		Matrix_Copy(U_minus, U_now, PARA_size);
		Matrix_Copy(U_now, U_plus, PARA_size);

		// ---------------------  collate data and exchange ---------------------------------------
    //MPI_Barrier(MPI_COMM_WORLD);
   	
   	if(size != 1){
      double Parcel[6] = {0};
      double post_box[6] = {0};

      // Populating the array that is shared by each process
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if(rank == 0){
          int counter = 0;
          for(int i = PARA_size-15; i < PARA_size-9; i++){
          	Parcel[counter] = U_now[i];
            counter++;
          } 
      }
      else if (rank == 1){
        int counter = 0;
        for(int i = 9; i < 15; i++){    // send 9 to 14
         	Parcel[counter] = U_now[i];
          counter++;
        }
      }

      // Both sending then receiving the data
      if(rank == 0){
        MPI_Send(Parcel, 6, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
        MPI_Recv(post_box, 6, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      else if(rank == 1){
        MPI_Recv(post_box, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(Parcel, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
      }

      // Process the data received in each processor
      if(rank == 0){
        int counter = 0;
        for(int i = PARA_size-6; i < PARA_size; i++){
          U_now[i] = post_box[counter]; 
          counter++;
        }
      }
      else if (rank == 1){
        int counter = 0;
        for(int i = 0; i < 6; i++){   // place 0 to 5
          U_now[i] = post_box[counter];
          counter++;
        }
      }
    }


    if(rank == 0 && (i < 50 || i > 9950) && output){
      std::cout << "\n-------------- Time step: " << i << " Tnow: " << t_now  << "-----------"; 
      std::cout << "\nU_now: " << std::endl;
      for(int i = 1; i < PARA_size; i+=3)
        std::cout << U_now[i] << "    ";
    }

		// -----------------------------------------------------------------

		delete[] result;
		delete[] Y2;
		delete[] Y3;
 		delete[] PARA_G_Fe;
 	}
 
  // ------------------- Task 4 (b) -------------------------------------
  /*
    Using process 1 to send its solution to process 0 to be written on the file,
    skipping the entire process if running off one process only
  */
  if(size != 1){
    MPI_Barrier(MPI_COMM_WORLD);
    double* Message = new double[PARA_size-15];
    if(rank == 1){
      int counter = 0;
      for(int i = 15; i < PARA_size; i++){
        Message[counter] = U_now[i];
        counter++;
      }
      MPI_Send(Message, PARA_size-15, MPI_DOUBLE, 0, 77, MPI_COMM_WORLD);    
    }
    if(rank == 0){ 
      double* Solution = new double[vec_size];
      initialise_dynamic_array(Solution, PARA_size);
      for(int i = 0; i < PARA_size; i++)
        Solution[i] = U_now[i];
      MPI_Recv(Message, PARA_size-15, MPI_DOUBLE, 1, 77, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      int counter = 0;
      for(int i = PARA_size; i < vec_size; i++){
        Solution[i] = Message[counter];
        counter++;
      }
      if(output)
        write_task4(Nx, Solution, vec_size, l);
    }
  }
  else{
    if(output){
      std::cout << "\nSolution for 1 processor only" << std::endl;
      for(int i = 1; i < vec_size; i+=3)
        std::cout << U_now[i] << "  ";
    }
  }
  if(output){
    std::cout << "FINAL SOL: " << U_now[int(Nx*n/2)] << std::endl;
  }

 
 	// use process 0 to write results to files

  // Clearing all the memory
  delete[] COEF1;
  delete[] COEF2;
  delete[] COEF3;

 	delete[] U_minus; 
 	delete[] U_now;
 	delete[] U_plus;

  delete[] PARA_G_Mm;
  delete[] PARA_G_Ke; 

	std::cout << "\n------------------- End of Task 4 ----------------------" << std::endl;
  return;
}
