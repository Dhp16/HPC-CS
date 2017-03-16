# include "task1.h"

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
void adding_three_arrays(double* result, double* A, double* B, double* C, const int Nx){
	int n = 3;
	for(int i = 0; i < Nx*n; i++){
		result[i] = A[i] + B[i] + C[i];
	}
}
void Diagonal_by_vector(double* result, const double* diagonal, const double* vector, const int Nx){
	const int n = 3;
	for(unsigned int i = 0; i < Nx*n; i++){
		result[i] = diagonal[i]*vector[i];
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
	int n = 3;
	int counter = 0;
	for(int i = k*Nx*n; i < (k+1)*Nx*n; i++){
		AA[i] = AA[i] + M[counter]; 
	}
	return;
}

void Build_Coef2(double* COEF2, const double G_Ke, const int Nx, const double delta_t){
	


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

void T2_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
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

	// COEF1: 1/delta_t^2 * [M]
	double* AA = new double[Nx*n];
	Global_Mass_Matrix(AA, Nx, rho, A, l);
	double u = 1/pow(delta_t, 2);
	Matrix_by_Scalar(AA, u, Nx);						// multiply Global Mass Matrix by the scalar
	Inverse_Diagonal_Matrix(AA, Nx);					// inverse the diagonal matrix
	double* AA2 = new double[Nx*n*9];
	Diagonal_to_Full_Size_Banded(AA2, AA, Nx);			// switch to banded form
	delete[] AA;
	double* COEF1 = new double[Nx*n*9];
	Matrix_Transformer(COEF1, AA2, Nx, k);				// switch to column major
	delete[] AA2;


	// COEF2: [K] - 2/dealtaT^2 [M]
	double* G_Me = new double[Nx*n];			
	Global_Mass_Matrix(G_Me, Nx, rho, A, l);			// generate M
	u = -2/pow(delta_t, 2);
	Matrix_by_Scalar(G_Me, u, Nx);						// multiply M by u 
	double* BB = new double[Nx*n*9];
	Global_Stiffness_Matrix(BB, Ke, Nx, k);				// generate & populate K
	Matrix_Add_Diagonal(BB, G_Me, Nx);					// add the two together
	double* COEF2 = new double[Nx*n*9];
	Matrix_Transformer(COEF2, BB, Nx, k);				// switch to column major
	delete[] BB;
														// at this stage COEF2 = [K] - 2[M]/delta_t^2  -    banded - column Major
	// COEF3: -1/delat_t^2[M]
	double* CC = new double[Nx*n];
	Global_Mass_Matrix(CC, Nx, rho, A, l);
	u = -1/pow(delta_t, 2);
	Matrix_by_Scalar(CC, u, Nx);						// COEF3 just diagonal
	double* CC2 = new double[Nx*n*Nx*n];
	Diagonal_to_Full_Size_Banded(CC2, CC, Nx);
	delete[] CC;
	double* COEF3 = new double[Nx*n*9];
	Matrix_Transformer(COEF3, CC2, Nx, k);
	delete[] CC2;

	// make size + transform
	
	// Initialising the integration scheme
	double* U_minus = new double[Nx*n];
	double* U_now = new double[Nx*n];
	double* U_plus = new double[Nx*n];
	initialise_dynamic_array(U_minus, Nx, n, 1);
	initialise_dynamic_array(U_now, Nx, n, 1);
	initialise_dynamic_array(U_plus, Nx, n, 1);


	//////////////////////////////////// testing ////////////////
	/*double t_now = 1*delta_t;
	std::cout << "delta_t: " << delta_t << std::endl;
	std::cout << "T: " << T << std::endl;
	double* Fi = new double[Nx*n];
	Global_Force_Vector(Fi, Fe, Nx, k);
	std::cout << "\nFi: " << std::endl;
	for(int i = 0; i < Nx*n; i++)
		std::cout << Fi[i] << " ";
	double K = t_now/T;
	std::cout << "\nK: " << K << std::endl;
	Matrix_by_Scalar(Fi, (double)t_now/T , Nx);*/

	////////////////////////////////////////////////////////////


	int t_steps = (int)std::round(Nt);					// need a fucntion to calculate the force every time step
	std::cout << "Starting loop" << std::endl;
	for(int i = 0; i < Nt; i++){
		std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
		// DGBMV
		// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
		double t_now = (i+1)*delta_t;
		double* Fi = new double[Nx*n];
		double* Y1 = new double[Nx*n];
		double* Y2 = new double[Nx*n];
		double* Y3 = new double[Nx*n];

		initialise_dynamic_array(Y1, Nx, n, 1);
		initialise_dynamic_array(Y2, Nx, n, 1);
		initialise_dynamic_array(Y3, Nx, n, 1);

		Global_Force_Vector(Fi, Fe, Nx, k, Fy);
		Matrix_by_Scalar(Fi, (double)t_now/T , Nx);

		int h = Nx*n;
		int kl = 4;
		int ku = 4;

		// COEF1*[Fi - COEF2*U_Now - COEF3*U_minus]
		std::cout << "\nFi " << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Fi[i] << " ";
		}

		std::cout << "\nU_now " << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << U_now[i] << " ";
		}

		F77NAME(dgbmv) ('n', h, h, kl, ku, 1.0, COEF2, h, U_now, 1, 0.0, Y2, 1);
		F77NAME(dgbmv) ('n', h, h, 1, 1, 1.0, COEF3, h, U_minus, 1, 0.0, Y3, 1);

		std::cout << "\nY2" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Y2[i] << " ";
		}
		std::cout << "\nY3" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Y3[i] << " ";
		}

		double* result = new double[Nx*n];
		adding_three_arrays(result, Fi, Y2, Y3, Nx);
		F77NAME(dgbmv) ('N', h, h, kl, ku, 1.0, COEF1, h, result, 1, 0.0, U_plus, 1);

		std::cout << "\nU_plus" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << U_plus[i] << " ";
		}		

		U_minus = U_now;
		U_now = U_plus;

		for(int i = 0; i < Nx*n; i++){
			U_minus[i] = U_now[i];
			U_now[i] = U_plus[i];
		}

		delete[] result;
		delete[] Fi;
		delete[] Y1;
		delete[] Y2;
		delete[] Y3;
	}
	// for(int i = 0; i < Nx*3; i++)

	delete[] COEF1;
	delete[] COEF2;
	delete[] COEF3;
	delete[] G_Me;
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

	// COEF1= 1/delta_t^2 * [M]
	double* M_inv = new double[Nx*n];
	Global_Mass_Matrix(M_inv, Nx, rho, A, l);
	double u = 1/pow(delta_t, 2);
	Matrix_by_Scalar(M_inv, u, Nx);						// multiply Global Mass Matrix by the scalar
	Inverse_Diagonal_Matrix(M_inv, Nx);					// inverse the diagonal matrix

	// COEF2= [K] - 2/dealtaT^2 [M]
	double* G_Me = new double[Nx*n];			
	Global_Mass_Matrix(G_Me, Nx, rho, A, l);			// generate M
	u = -2/pow(delta_t, 2);
	Matrix_by_Scalar(G_Me, u, Nx);						// multiply M by u 
	double* BB = new double[Nx*n*9];
	Global_Stiffness_Matrix(BB, Ke, Nx, k);				// generate & populate K
	Matrix_Add_Diagonal(BB, G_Me, Nx);					// add the two together
	double* COEF2 = new double[Nx*n*9];
	Matrix_Transformer(COEF2, BB, Nx, k);	
	
	// switch to column major
	//delete[] BB;
														// at this stage COEF2 = [K] - 2[M]/delta_t^2  -    banded - column Major
	// COEF3= -1/delat_t^2[M]
	double* COEF3 = new double[Nx*n];
	Global_Mass_Matrix(COEF3, Nx, rho, A, l);
	u = -1/pow(delta_t, 2);
	Matrix_by_Scalar(COEF3, u, Nx);						// COEF3 just diagonal
	
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
		std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
		// DGBMV
		// http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
		double t_now = (i+1)*delta_t;
		double* Fi = new double[Nx*n];
		double* Y1 = new double[Nx*n];
		double* Y2 = new double[Nx*n];
		double* Y3 = new double[Nx*n];

		initialise_dynamic_array(Y1, Nx, n, 1);
		initialise_dynamic_array(Y2, Nx, n, 1);
		initialise_dynamic_array(Y3, Nx, n, 1);
		Global_Force_Vector(Fi, Fe, Nx, k, Fy, t_now/T); 	// populating the force vector and scaling the forces by t_now/T

		int h = Nx*n;
		int kl = 4;
		int ku = 4;

		// COEF1*[Fi - COEF2*U_Now - COEF3*U_minus]
		std::cout << "\nFi " << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Fi[i] << " ";
		}

		std::cout << "\nU_now " << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << U_now[i] << " ";
		}

		F77NAME(dgbmv) ('n', h, h, kl, ku, 1.0, COEF2, h, U_now, 1, 0.0, Y2, 1);
		F77NAME(dgbmv) ('n', h, h, 0, 0, 1.0, COEF3, 1, U_minus, 1, 0.0, Y3, 1);

		std::cout << "\nY2" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Y2[i] << " ";
		}
		std::cout << "\nY3" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << Y3[i] << " ";
		}

		double* result = new double[Nx*n];
		adding_three_arrays(result, Fi, Y2, Y3, Nx);

		int lda = 1;
		F77NAME(dgbmv) ('n', h, h, 0, 0, 1.0, M_inv, 1, result, 1, 0.0, U_plus, 1);

		std::cout << "\nU_plus" << std::endl;
		for(int i = 0; i < Nx*n; i++){
			std::cout << U_plus[i] << " ";
		}		

		U_minus = U_now;
		U_now = U_plus;

		delete[] result;
		delete[] Fi;
		delete[] Y1;
		delete[] Y2;
		delete[] Y3;
	}

	delete[] BB;
	delete[] M_inv;
	delete[] COEF2;
	delete[] COEF3;
	delete[] G_Me;
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
	
	std::cout << "COEF2" << std::endl;
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
	Matrix_Transformer(COEF2, BB, Nx);					// Wrong 
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




























