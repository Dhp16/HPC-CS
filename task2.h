/* 
 Imperial College London
 HPC Assignment Task 3
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 2

 The variables used are all explained in function T2_inputs

 DGBMV documentation
 http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
*/

// ---------------------------------- Task 2 (b) ------------------
// Build Global Mass Matrix
void Global_Mass_Matrix(double* G_Mm, const int Nx, const double rho, const double A, const double l){					 
	// Populates the G_Mm (Global Mass Matrix) node by node according to Equation 3
	const int n = 3;								
	const double alpha = 1./24; 
	for(int i = 0; i < Nx*n; i+=3){
		G_Mm[i] = 1.;
		G_Mm[i+1] = 1.;
		G_Mm[i+2] = 2*alpha*pow(l, 2);
	}
	for(int i = 0; i < Nx*n; i++){					// multiply each element by rho * A * l
		G_Mm[i] = G_Mm[i]*rho*A*l;
	}
	return;
}
// -------------------------------- Task 2 (c) -------------------
// COEF1 = delta_t^2/[M]
void Build_Coef1(double* COEF1, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF1[i] = pow(delta_t, 2)/G_Mm[i];
	return;
}
// COEF2 = [K] - 2/delta_t^2 *[M]
void Build_Coef2(double* COEF2, const double* G_Ke, const double* G_Mm, const int Nx, const double delta_t){
	int n = 3;
	for(int i = 0; i < Nx*n*9; i++)
		COEF2[i] = G_Ke[i];
	int counter = 0;
	for(int i = 4*Nx*n; i < 5*Nx*n; i++){
		COEF2[i] = (COEF2[i] - 2/(pow(delta_t, 2))*G_Mm[counter]);
		counter++;
	}
	Matrix_Transformer_Imp(COEF2, Nx, 0);
}
// COEF3 = -1/deltat^2 * [M]
void Build_Coef3(double* COEF3, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF3[i] = -G_Mm[i]/pow(delta_t, 2);
	return; 	
}
// ------------------------------- Task 2 (e) -------------------- 
// Writting mid point deflection @ each time step
void write_task2_displacements(std::vector<double> tsteps, std::vector<double> dispMid){	
	std::ofstream file;
	std::string file_name = "Task2_Mid_point_vs_Time.txt"; 
	file.open(file_name);
	for(int i = 0; i < dispMid.size(); i++)
		file << tsteps[i] << std::endl; 
	for(int i = 0; i < dispMid.size(); i++)
		file << dispMid[i] << std::endl; 
	file.close();
	return;
}
// ------------------------------- Task 2 (f) --------------------
// Writting amplitude of mid-point oscillation for each loading time
void write_oscillations(const std::vector<double> amplitudes, const std::vector<double> lt){
	if(amplitudes.size() != lt.size()){
		std::cout << "ERROR: Vectors of unequal size, cancelling write to file for oscillations." << std::endl;
		return;
	}
	std::ofstream file;
	std::string file_name = "Task3_Oscillations.txt"; 
	file.open(file_name);
	for(int i = 0; i < amplitudes.size(); i++)
		file << std::setprecision(12) << amplitudes[i] << std::endl; 
	for(int i = 0; i < lt.size(); i++)
		file << std::setprecision(12) << lt[i] << std::endl; 
	file.close();
	return;
}

// -------------------------- Primary Function -------------------
void T2_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){ 

	/*
		The explicit scheme was solved using inversion of the diagonal matrix [M] rather than using a lapack 

	*/

	std::cout << "\n---------------------- Task 2 ------------------------" << std::endl;
	const int n = 3;
	int vec_size = Nx*n; 
	double delta_t = T/Nt;
	int sizeA = 10;
 
	double Fe[6] = {0};
	pp_Fe(Fe, l, Nx, Qy, Qx);						// populate elemental force vector
	int k = 0;
	double Ke[36] = {0};	
	pp_Ke(Ke, A, E, I, l);							// populate elemental stiffness matrix

	double* G_Mm = new double[vec_size];			// Declare the Global Mass Matrix
	Global_Mass_Matrix(G_Mm, Nx, rho, A, l);		// Populate the Global Mass Matrix

	double* G_Ke = new double[vec_size*9];			// Declare the Global Stiffness Matrix
	Global_Stiffness_Matrix(G_Ke, Ke, Nx, k);		// Populate the Global Stiffness Matrix

	double* COEF1 = new double[vec_size];			
	double* COEF2 = new double[vec_size*9];			
	double* COEF3 = new double[vec_size];			

	// COEF1= 1/delta_t^2 * [M]
	Build_Coef1(COEF1, G_Mm, Nx, delta_t);			// COEF1 = delta_t^2/[M]
	Build_Coef2(COEF2, G_Ke, G_Mm, Nx, delta_t);    // COEF2 = [K] - 2/delta_t^2 *[M]
	Build_Coef3(COEF3, G_Mm, Nx, delta_t);			// COEF3 = -1/deltat^2 * [M] 

	// Control of loading time
	/* 
		As it stands the load will be applied linearly from [0; T], however, changing the "factor" below to a higher number will decrease the loading time. A vector of 4000 different 
		factors were used to plot the results shown for amplitude against loading time in the report.
	*/
	double factor = 1.0;


	// Initialising the integration scheme
	double* U_minus = new double[Nx*n];
	double* U_now = new double[Nx*n];
	double* U_plus = new double[Nx*n];

	initialise_dynamic_array(U_minus, vec_size);
	initialise_dynamic_array(U_now, vec_size);
	initialise_dynamic_array(U_plus, vec_size); 	
   
	std::vector<double> tsteps;
	std::vector<double> dispMid;   
  
	bool loading_finished = false;
	double loading_time = 0;
	double end_loading_time_step = 0;
	int t_steps = (int)std::round(Nt);								// need a fucntion to calculate the force every time step

	for(int i = 0; i < Nt; i++){
		double t_now = i*delta_t; 
		tsteps.push_back(t_now); 
		double* Fi = new double[Nx*n]; 
		double* Y2 = new double[Nx*n];
		double* Y3 = new double[Nx*n]; 
		initialise_dynamic_array(Y2, vec_size);
		initialise_dynamic_array(Y3, vec_size);
		initialise_dynamic_array(Fi, vec_size);
		if(factor*t_now/T >= 1.0)
			Global_Force_Vector(Fi, Fe, Nx, k, Fy, 1); 					// populating the force vector without any scaling as the maximum load has been reached
 		else 
 			Global_Force_Vector(Fi, Fe, Nx, k, Fy, factor*t_now/T);		// populating the force vector and scaling the forces by t_now/T
	
 		if(factor*t_now/T >= 1.0 && loading_finished == false){			// flag the fact the maximum loading has been reached
 			loading_finished = true;
 			loading_time = t_now;  
 			end_loading_time_step = i; 
 		} 

		double* copy_U_now = new double[Nx*n]; 
		Matrix_Copy(copy_U_now, U_now, vec_size);														// Solving the System
		F77NAME(dgbmv)('n', vec_size, vec_size, 4, 4, 1.0, COEF2, 9, copy_U_now, 1, 0.0, Y2, 1);		// Y2 = COEF2*U_n
		delete[] copy_U_now;
		Diagonal_by_Vector(Y3, COEF3, U_minus, vec_size);												// Y3 = COEF3*U_n-1
		double* result = new double[Nx*n];
		initialise_dynamic_array(result, Nx*n);
		adding_three_arrays(result, Fi, Y2, Y3, vec_size);

		Diagonal_by_Vector(U_plus, COEF1, result, vec_size);	

		dispMid.push_back(U_plus[(Nx*3)/2]);															// store the displacement of the mid point for each time step (task 2 (e))

		Matrix_Copy(U_minus, U_now, vec_size);															// Switching the matrices for the next time step
		Matrix_Copy(U_now, U_plus, vec_size);
		delete[] result; 
		delete[] Fi;
		delete[] Y2; 
		delete[] Y3; 

		if(i != Nt-1) 
			initialise_dynamic_array(U_plus, Nx*3);
	}
	
	write_task_solution(U_now, Nx*n, "2");																// Writting Proof of Solution for Task 2 (c) 
	//write_task2_displacements(tsteps, dispMid);														// Writting the data needed for Task 2 (e)

	// Clearing memory
	delete[] U_minus;
	delete[] U_plus;
	delete[] U_now;

	delete[] COEF1;
	delete[] COEF2;
	delete[] COEF3;
	delete[] G_Mm;
	std::cout << "------------------- End of Task 2 ----------------------" << std::endl;
	return;
}

/* 
	As it stands the load will be applied linearly from [0; T], however, changing the "factor" below to a higher number will decrease the loading time. A vector of 4000 different 
	factors were used to plot the results shown for amplitude against loading time in the report.
*/








