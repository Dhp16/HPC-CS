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
	const int n = 3;
	const double alpha = 1./24; 
	for(int i = 0; i < Nx*n; i+=3){
		G_Mm[i] = 1.;
		G_Mm[i+1] = 1.;
		G_Mm[i+2] = 2*alpha*pow(l, 2);
	}
	for(int i = 0; i < Nx*n; i++){
		G_Mm[i] = G_Mm[i]*rho*A*l;
	}
	return;
}
// -------------------------------- Task 2 (c) -------------------
void Build_Coef1(double* COEF1, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF1[i] = pow(delta_t, 2)/G_Mm[i];
	return;
}
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
void Build_Coef3(double* COEF3, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF3[i] = -G_Mm[i]/pow(delta_t, 2);
	return; 	
}
// ------------------------------- Task 2 (e) --------------------
void write_task2(std::string file_name, std::vector<double> tsteps, std::vector<double> dispMid){
	std::ofstream file;
	file_name = "Task2_" + file_name + ".txt"; 
	file.open(file_name);
	for(int i = 0; i < dispMid.size(); i++)
		file << tsteps[i] << std::endl; 
	for(int i = 0; i < dispMid.size(); i++)
		file << dispMid[i] << std::endl; 
	file.close();
	return;
}
// ------------------------------- Task 2 (f) --------------------
void write_oscillations(const std::vector<double> amplitudes, const std::vector<double> lt){
	if(amplitudes.size() != lt.size()){
		std::cout << "ERROR: Vectors of unequal size, cancelling write to file for oscillations." << std::endl;
		return;
	}
	std::ofstream file;
	std::string file_name = "Task2_Oscillations.txt"; 
	file.open(file_name);
	for(int i = 0; i < amplitudes.size(); i++)
		file << amplitudes[i] << std::endl; 
	for(int i = 0; i < lt.size(); i++)
		file << lt[i] << std::endl; 
	file.close();
	return;
}


void T2_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){ 

	/*
		The explicit scheme was solved using inversion of the diagonal matrix [M] rather than using a lapack 

	*/


	std::cout << "---------------------- Task 2 ------------------------" << std::endl;
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

	double* COEF1 = new double[vec_size];			// COEF1 = delta_t^2/[M]
	double* COEF2 = new double[vec_size*9];			// COEF2 = [K] - 2/delta_t^2 *[M]
	double* COEF3 = new double[vec_size];			// COEF3 = -1/deltat^2 * [M] 

	// COEF1= 1/delta_t^2 * [M]
	Build_Coef1(COEF1, G_Mm, Nx, delta_t);
	Build_Coef2(COEF2, G_Ke, G_Mm, Nx, delta_t);       
	Build_Coef3(COEF3, G_Mm, Nx, delta_t);

	// Initialising the integration scheme
	double* U_minus = new double[Nx*n];
	double* U_now = new double[Nx*n];
	double* U_plus = new double[Nx*n];  
	initialise_dynamic_array(U_minus, Nx, n, 1);
	initialise_dynamic_array(U_now, Nx, n, 1);
	initialise_dynamic_array(U_plus, Nx, n, 1); 	
 
	std::vector<double> tsteps;
	std::vector<double> dispMid;  
 
	std::vector<double> factors;
	double fac = 1.001;		
	for(int k = 0; k < 5000; k++){
		factors.push_back(fac);
		fac += 0.002;
	}

	std::vector<double> amplitudes;
	std::vector<double> lt;

	std::cout << "starting loop" << std::endl;
	for(int k = 0; k < factors.size(); k++){
		if(k % 100 == 0)
			std::cout << k << std::endl;
		double factor = factors[k];						// time factor to control loading time

		bool loading_finished = false;
		int t_steps = (int)std::round(Nt);					// need a fucntion to calculate the force every time step
		//std::cout << "\nStarting loop" << std::endl;
		for(int i = 0; i < Nt; i++){
			double t_now = i*delta_t;
			tsteps.push_back(t_now);
			double* Fi = new double[Nx*n];
			double* Y2 = new double[Nx*n];
			double* Y3 = new double[Nx*n]; 
			initialise_dynamic_array(Y2, Nx, n, 1);
			initialise_dynamic_array(Y3, Nx, n, 1);
			initialise_dynamic_array(Fi, Nx, n, 1);
			if(factor*t_now/T > 1)
				Global_Force_Vector(Fi, Fe, Nx, k, Fy, 1); 	// populating the force vector and scaling the forces by t_now/T
 			else 
 				Global_Force_Vector(Fi, Fe, Nx, k, Fy, factor*t_now/T);
	
 			if(factor*t_now/T >= 1 && loading_finished == false){
 				//std::cout << "End of loading: " << t_now << std::endl;
 				//lt.push_back(t_now); 
 				loading_finished = true;
 			} 
																						// COEF1*[Fi - COEF2*U_Now - COEF3*U_minus]

			double* copy_U_now = new double[Nx*n];
			Matrix_Copy(copy_U_now, U_now, vec_size);

			F77NAME(dgbmv)('n', vec_size, vec_size, 4, 4, 1.0, COEF2, 9, copy_U_now, 1, 0.0, Y2, 1);		// Y2 = COEF2*U_n
			delete[] copy_U_now;

			Diagonal_by_Vector(Y3, COEF3, U_minus, vec_size);										// Y3 = COEF3*U_n-1

			double* result = new double[Nx*n];
			initialise_dynamic_array(result, Nx, n, 1);
			adding_three_arrays(result, Fi, Y2, Y3, vec_size);
	
			Diagonal_by_Vector(U_plus, COEF1, result, vec_size);

			dispMid.push_back(U_plus[((Nx/2)*3) - 2]);

			Matrix_Copy(U_minus, U_now, vec_size);
			Matrix_Copy(U_now, U_plus, vec_size);

			delete[] result; 
			delete[] Fi;
			delete[] Y2; 
			delete[] Y3;

			if(i != Nt-1) 
				initialise_dynamic_array(U_plus, Nx, n, 1);
		}

		int DMsize = (dispMid.size()/factor)+1;

		double max = -100;
		double min = 100;
		for(int i = DMsize; i < dispMid.size(); i++){
			if(dispMid[i] > max)
				max = dispMid[i];
			if(dispMid[i] < min)
				min  = dispMid[i];
		}

		double amplitude = max - min;
		//std::cout << "amplitude: " << amplitude << std::endl;
		amplitudes.push_back(amplitude);
		//std::cout << "max: " << max << "  min: " << min << "  amplitude: " << amplitude << std::endl; 
		bool full_time = false;

		//std::string file_name;
		//if(full_time)
		//	file_name = "FT";
		//else 
		//	file_name = "HT";
	}

	//std::cout << "lt size: " << lt.size() << std::endl;
	//std::cout << "amplitudes: " << amplitudes.size() << std::endl;
	write_oscillations(amplitudes, factors);

	//write_task2(file_name, tsteps, dispMid);

	// cannot delete U_Plus & U_Minus


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










