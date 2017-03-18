// DGBMV
// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888



// ---------------------------------- Tools for Task 2 ------------
void Diagonal_to_Full_Size_Banded(double* output, const double* input, const int Nx){
	int n = 3;
	int counter = 0;
	for(int i = 0; i < Nx*n*9; i++){
		output[i] = 0;
	}
	for(int i = 4*Nx*n; i < 5*Nx*n; i++){
		output[i] = input[counter];
		counter++;
	}
	return;
}
void Diagonal_by_vector_Imp(double* diagonal, const double* vector, const int Nx){
	const int n = 3;
	for(unsigned int i = 0; i < Nx*n; i++){
		diagonal[i] = diagonal[i]*vector[i];
	}
	return;
}

// ---------------------------------- Task 2 (b) ------------------
// Build Global Mass Matrix
void Global_Mass_Matrix(double* G_Mm, const int Nx, const double rho, const double A, const double l){					
	// creating just the diagonal 
	int n = 3;
	double alpha = 1./24; 
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
void Matrix_by_Scalar(double* M, const double u, const int Nx){
	const int n = 3;
	for(int i = 0; i < Nx*n; i++){
		M[i] = M[i]*u;
	}
	return;
}
void Matrix_Add_Diagonal(double* AA, const double* M, const int Nx, const int k = 4){
	// k: number of upper diagonals
	int n = 3;
	int counter = 0;
	for(int i = k*Nx*n; i < (k+1)*Nx*n; i++){
		AA[i] = AA[i] + M[counter]; 
	}
	return;
}

void Build_Coef1(double* COEF1, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF1[i] = pow(delta_t, 2)/G_Mm[i];
	return;
}
void Build_Coef2(double* COEF2, const double* G_Ke, const double* G_Mm, const int Nx, const double delta_t){
	int n = 3;
	for(int i = 0; i < Nx*n*9; i++)
		COEF2[i] = -G_Ke[i];
	int counter = 0;
	for(int i = 4*Nx*n; i < 5*Nx*n; i++){
		COEF2[i] = COEF2[i] + 2/pow(delta_t, 2)*G_Mm[counter];
		counter++;
	}
	Matrix_Transformer_Imp(COEF2, Nx, 0);
}
void Build_Coef3(double* COEF3, const double* G_Mm, const int Nx, const double delta_t){
	for(int i = 0; i < Nx*3; i++)
		COEF3[i] = -G_Mm[i]/pow(delta_t, 2);
	return; 	
}
 
// Explicit time integration (Central Difference Method):
void Explicit_time_integration(double* G_Mm, const int Nx, double){
}
void Inverse_Diagonal_Matrix(double* M, const int Nx){
	// M any diagonal matrix
	// inverse of a diagonal matrix is obtained by replacing each element by its inverse
	int n = 3;
	for(int i = 0; i < Nx*n; i++){
		M[i] = 1./M[i];
	}
	return;
}

void T2_inputsV2(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){ 

	std::cout << "\nStarting Task 2" << std::endl;

	const int n = 3; 
	double delta_t = T/Nt;
	int sizeA = 10;
 
	double Fe[6] = {0};
	pp_Fe(Fe, l, Nx, Qy, Qx);
	int k = 0;
	double Ke[36] = {0};
	pp_Ke(Ke, A, E, I, l);

	double* G_Mm = new double[Nx*n];
	Global_Mass_Matrix(G_Mm, Nx, rho, A, l);
	double* G_Ke = new double[Nx*n*9];
	Global_Stiffness_Matrix(G_Ke, Ke, Nx, k);	// check this

	double* COEF1 = new double[Nx*n];			// COEF1 = delta_T^2/[M]
	double* COEF2 = new double[Nx*n*9];			// COEF2 = [K] - 2/dealtaT^2 [M]
	double* COEF3 = new double[Nx*n];			// COEF3 = -1/deltat^2 * [M] 

	// COEF1= 1/delta_t^2 * [M]
	Build_Coef1(COEF1, G_Mm, Nx, delta_t);
	Build_Coef2(COEF2, G_Ke, G_Mm, Nx, delta_t);        // at this stage COEF2 = [K] - 2[M]/delta_t^2  -    banded - column Major
	Build_Coef3(COEF3, G_Mm, Nx, delta_t);


	std::cout << "\nCOEF1" << std::endl;
	for(int i = 0; i < Nx*n; i++){
		std::cout << std::setw(10) << COEF1[i] << "         ";
	}

	std::cout << "\nCOEF3" << std::endl;
	for(int i = 0; i < Nx*n; i++){
		std::cout << std::setw(10) << COEF3[i] << "         ";
	}

	std::cout << "\nCOEF2" << std::endl;
	for(int i = 0; i < 9*Nx*n; i++){
		if(i % 9 ==0)
			std::cout << "\n";
		std::cout << std::setw(10) << COEF2[i] << "         ";
	}

	// Initialising the integration scheme
	double* U_minus = new double[Nx*n];
	double* U_now = new double[Nx*n];
	double* U_plus = new double[Nx*n];
	initialise_dynamic_array(U_minus, Nx, n, 1);
	initialise_dynamic_array(U_now, Nx, n, 1);
	initialise_dynamic_array(U_plus, Nx, n, 1);

	int t_steps = (int)std::round(Nt);					// need a fucntion to calculate the force every time step
	std::cout << "\nStarting loop" << std::endl;
	for(int i = 0; i < Nt; i++){
		// DGBMV
		// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
		double t_now = i*delta_t;
		double* Fi = new double[Nx*n];
		double* Y1 = new double[Nx*n];
		double* Y2 = new double[Nx*n];
		double* Y3 = new double[Nx*n];
		initialise_dynamic_array(Y1, Nx, n, 1);
		initialise_dynamic_array(Y2, Nx, n, 1);
		initialise_dynamic_array(Y3, Nx, n, 1);
		initialise_dynamic_array(Fi, Nx, n, 1);
		Global_Force_Vector(Fi, Fe, Nx, k, Fy, t_now/T); 	// populating the force vector and scaling the forces by t_now/T
 
		int h = Nx*n;
		int kl = 4;
		int ku = 4;

		// COEF1*[Fi - COEF2*U_Now - COEF3*U_minus]

		double* copy_U_now = new double[Nx*n];
		for(int i = 0; i < Nx*n; i++)
			copy_U_now[i] = U_now[i];

		F77NAME(dgbmv)('n', h, h, kl, ku, 1.0, COEF2, h, copy_U_now, 1, 0.0, Y2, 1);		// Y2 = COEF2*U_n
		delete[] copy_U_now;

		Diagonal_by_vector(Y3, COEF3, U_minus, Nx);										// Y3 = COEF3*U_n-1

		double* result = new double[Nx*n];
		initialise_dynamic_array(result, Nx, n, 1);
		adding_three_arrays(result, Fi, Y2, Y3, Nx);

		Diagonal_by_vector(U_plus, COEF1, result, Nx);

		for(int i = 0; i < Nx*n; i++){
			U_minus[i] = U_now[i];
			U_now[i] = U_plus[i];
		}	 

		if(i < 10){
			std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
			std::cout << "\nY2" << std::endl;
			for(int i = 0; i < Nx*n; i++){
				std::cout << Y2[i] << " ";
			}
			std::cout << "\nY3" << std::endl;
			for(int i = 0; i < Nx*n; i++){
				std::cout << Y3[i] << " ";
			}
			std::cout << "\nU_plus" << std::endl;
			for(int i = 1; i < Nx*n; i+=3){
				std::cout << U_plus[i] << " ";
			}		 
		}

		delete[] result;
		delete[] Fi;
		delete[] Y1;
		delete[] Y2;
		delete[] Y3;
	}

	std::cout << "\n\n\n------ Final Displacement ------" << std::endl;
	for(int i = 1; i < Nx*n; i+=3){
		std::cout << U_plus[i] << " ";
	}		

	// cannot delete U_Plus & U_Minus

	delete[] U_minus;
	delete[] U_plus;
	delete[] U_now;

	delete[] COEF1;
	delete[] COEF2;
	delete[] COEF3;
	delete[] G_Mm;
	return;
}

void T2_inputsV3(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){ 

	std::cout << "\nStarting Task 2" << std::endl;
	const int n = 3;
	double delta_t = T/Nt;
	int sizeA = 10;

	double Fe[6] = {0};
	pp_Fe(Fe, l, Nx, Qy, Qx);
	int k = 0;
	double Ke[36] = {0};
	pp_Ke(Ke, A, E, I, l);

	// COEF1= 1/delta_t^2 * [M]
	double* M_inv = new double[Nx*n];
	Global_Mass_Matrix(M_inv, Nx, rho, A, l);
	double u = 1/pow(delta_t, 2);
	Matrix_by_Scalar(M_inv, u, Nx);						// multiply Global Mass Matrix by the scalar
	Inverse_Diagonal_Matrix(M_inv, Nx);					// inverse the diagonal matrix

	double* U_i = new double[Nx*n];
	initialise_dynamic_array(U_i, Nx, n, 1);
	double* U_i_minus = new double[Nx*n];
	initialise_dynamic_array(U_i_minus, Nx, n, 1);
	double* U_i_plus = new double[Nx*n];
	initialise_dynamic_array(U_i_plus, Nx, n, 1);
	
	// ------------------  COEF2= [K] - 2/dealtaT^2 [M] ------------------------
	double* G_Me = new double[Nx*n];			
	Global_Mass_Matrix(G_Me, Nx, rho, A, l);			// generate M
	u = -2/pow(delta_t, 2);
	Matrix_by_Scalar(G_Me, u, Nx);						// multiply M by u 
	double* BB = new double[Nx*n*9];
	Global_Stiffness_Matrix(BB, Ke, Nx);				// generate & populate K
	Matrix_Add_Diagonal(BB, G_Me, Nx, 4);				// add the two together
	delete G_Me;
	double* COEF2 = new double[Nx*n*9];
	output_array(BB, "Coef2", Nx*3, 9);
	Matrix_Transformer(COEF2, BB, Nx);					//// Wrong 
	delete[] BB;
	// -------------------------------------------------------------------------
 
	// ------------------ COEF3= -1/2delta_t^2 [M] -----------------------------
	double* COEF3 = new double[Nx*n];
	Global_Mass_Matrix(COEF3, Nx, rho, A, l);			// correcting previous changes to G_Me
	u = -1/pow(delta_t, 2);
	Matrix_by_Scalar(COEF3, u, Nx);		
	// -------------------------------------------------------------------------

	int t_steps = (int)std::round(Nt);					// need a fucntion to calculate the force every time step
	std::cout << "\nStarting loop" << std::endl;
	for(int i = 1; i < Nt+1; i++){
		
		// declaring variables to store the results of calculations
		double* Y2 = new double[Nx*n];		// COEF2*u{i}
		double* Y3 = new double[Nx*n];		// COEF2*u{i-1}
		double* RHS = new double[Nx*n];		

		// Calculate F
		double t_now = i*delta_t; 
		double* Fi = new double[Nx*n];
		Global_Force_Vector(Fi, Fe, Nx, k, Fy, t_now/T); 	// populating the force vector and scaling the forces by t_now/T

		int lda = 9; // number of rows
		F77NAME(dgbmv) ('n', Nx*n, Nx*n, 4, 4, 1.0, COEF2, lda, U_i, 1, 0.0, Y2, 1);
		//F77NAME(dgbmv) ('n', Nx*n, Nx*n, 0, 0, 1.0, COEF3, 1, U_i_minus, 1, 0.0, Y3, 1);
		Diagonal_by_vector(Y3, COEF3, U_i_minus, Nx);

		adding_three_arrays(RHS, Fi, Y2, Y3, Nx);

		// mutiplying by the inverse:
		//F77NAME(dgbmv) ('n', Nx*n, Nx*n, 0, 0, 1.0, M_inv, Nx*n, RHS, 1, 0.0, U_i_plus, 1);
		Diagonal_by_vector(U_i_plus, M_inv, RHS, Nx);

		if(i < 5){
			std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
			std::cout << "\nY2: " << std::endl;
			for(int i = 0; i < Nx*n; i++){
				std::cout << " " <<Y2[i] << std::endl;
			}
			std::cout << "\nY3: " << std::endl;
			for(int i = 0; i < Nx*n; i++){
				std::cout << " " <<Y3[i] << std::endl;
			}
			std::cout << "\nRHS: " << std::endl; 
			for(int i = 0; i < Nx*n; i++){
				std::cout << " " <<RHS[i] << std::endl;
			}
			std::cout << "\nDisplacement" << std::endl;
			for(int i = 1; i < Nx*n; i+=3){
				std::cout << " " << U_i[i] << std::endl;
			}

		}

		// assigning the values for the following time step
		for(int i = 0; i < Nx*n; i++){
			U_i_minus[i] = U_i[i];
			U_i[i] = U_i_plus[i];
		}

		delete[] Fi;
		delete[] Y2;
		delete[] Y3;
		delete[] RHS;
	} 
	std::cout << "\nFinal Displacement at ts: "<< Nt << std::endl;
	for(int i = 1; i < Nx*n; i+=3)
			std::cout << " " << U_i[i] << std::endl;
	

	delete[] M_inv;
	delete[] COEF2;
	delete[] COEF3;
	delete[] U_i;
	delete[] U_i_plus;
	delete[] U_i_minus;
	return;
}




























