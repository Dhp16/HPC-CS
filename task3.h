


void K_Effective(double* K_effective, const double* G_Mm, const int Nx, const double delta_t){
	// assuming the K_effective that is input is already equal to G_Ke
	int n = 3;
	double* G_Mm_copy = new double[Nx*n];
	Matrix_Copy(G_Mm_copy, G_Mm, Nx*n);
	double u = 1/(0.250*delta_t);
	Matrix_by_Scalar(G_Mm_copy, u, Nx);
	Matrix_Add_Diagonal(K_effective, G_Mm_copy, Nx, 0);
	delete[] G_Mm_copy;
}

void T3_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){

	double beta = 0.250; 
	double gamma = 0.5;

	int n = 3;
	int vec_size = Nx*n;
	double delta_t = T/Nt;											// time step

	double Fe[6] = {0};											// elemental force vector
	pp_Fe(Fe, l, Nx, Qy, Qx);									// populate elemental force vector
	int k = 0;
	double Ke[36] = {0};				
	pp_Ke(Ke, A, E, I, l);										// populate elemental stiffness matrix

	double* G_Mm = new double[vec_size];
	Global_Mass_Matrix(G_Mm, Nx, rho, A, l);
	double* G_Ke = new double[vec_size*9];
	Global_Stiffness_Matrix(G_Ke, Ke, Nx, k);

	double* K_effective = new double[Nx*n*13];
	Global_Stiffness_Matrix(K_effective, Ke, Nx, 4);
	K_Effective(K_effective, G_Mm, Nx, delta_t);		// inside the loop beacuse DGBSV mangles with it

	double* DIS = new double[Nx*n];
	double* VEL = new double[Nx*n];
	double* ACC = new double[Nx*n];
	initialise_dynamic_array(DIS, Nx, n, 1);
	initialise_dynamic_array(VEL, Nx, n, 1);
	initialise_dynamic_array(ACC, Nx, n, 1);
	double* DIS_plus = new double[Nx*n];
	double* VEL_plus = new double[Nx*n];
	double* ACC_plus = new double[Nx*n];
	initialise_dynamic_array(DIS_plus, Nx, n, 1);
	initialise_dynamic_array(VEL_plus, Nx, n, 1);
	initialise_dynamic_array(ACC_plus, Nx, n, 1);

	for(int i = 0; i < Nt; i++){
		double t_now = i*delta_t;

		double* Fi = new double[Nx*n];
		initialise_dynamic_array(Fi, Nx, n, t_now/T);
		Global_Force_Vector(Fi, Fe, Nx, k, Fy, 1);

		// Equation 12
		double* DIS12 = new double[vec_size]; Matrix_Copy(DIS12, DIS, vec_size);
		double* VEL12 = new double[vec_size]; Matrix_Copy(VEL12, VEL, vec_size);
		double* ACC12 = new double[vec_size]; Matrix_Copy(ACC12, ACC, vec_size);
		
		vector_by_scalar(DIS12, 1/(beta*pow(delta_t, 2)), vec_size);					// willl need to use some of these again without the scalling
		vector_by_scalar(VEL12, 1/(beta*delta_t), vec_size);
		vector_by_scalar(ACC12, 1/(2*beta)-1, vec_size);

		double* A = new double[vec_size];
		adding_three_arrays(A, DIS12, VEL12, ACC12, vec_size);
		double* B = new double[vec_size];
		Diagonal_by_Vector(B, G_Mm, A, vec_size);
		delete[] A;
		adding_to_array(B, Fi, vec_size);															// B = {F}n+1  + [M]{......}
		int* ipiv = new int[Nx*3];
    	int info = 1;

    	double* K_effective_copy = new double[vec_size*13];
    	Matrix_Copy(K_effective_copy, K_effective, vec_size*13);								// protecting K_effective from the changes incureed by dgbsv
    	Matrix_Transformer_Imp(K_effective_copy, Nx, 4);


    	F77NAME(dgbsv) (vec_size, 4, 4, 1, K_effective, 1+2*4+4, ipiv, B, vec_size, &info);	// result should be ported to DIS_plus10--> save KEff for next loop
    	if(i == 0 && info != 0)
    		std::cout << "dgbsv failed!" << std::endl;
    	delete[] K_effective_copy;
    	Matrix_Copy(DIS_plus, B, vec_size);
    	delete[] B;

    	delete[] DIS12;
    	delete[] VEL12;
    	delete[] ACC12;

		// Equation 10:
    	double* DIS10 = new double[Nx*n]; 		Matrix_Copy(DIS10, DIS, vec_size);
    	double* DIS_plus10 = new double[Nx*n];  Matrix_Copy(DIS_plus10, DIS_plus, vec_size);
		double* VEL10 = new double[Nx*n]; 		Matrix_Copy(VEL10, VEL, vec_size);
		double* ACC10 = new double[Nx*n]; 		Matrix_Copy(ACC10, ACC, vec_size);

    	sub_from_vector(DIS_plus10, DIS, vec_size);											// DIS_plus = DIS_plus - DIS_plus--> {u}n+1 = ({u}n+1 - {u}n)
    	vector_by_scalar(DIS_plus10, 1/(beta*pow(delta_t, 2)), vec_size);				// DIS_plus = DIS_plus*1/(beta*delta_t^2)		

    	vector_by_scalar(VEL10, -1/(beta*delta_t), vec_size);
    	vector_by_scalar(ACC10, 1-1/(2*beta), vec_size);
    	adding_three_arrays(ACC_plus, DIS_plus10, VEL10, ACC10, vec_size);	// ACC_plus = 1/(beta*delta_t^2)*({u}n+1 - {u}n) - 1/(beta*delta_t)*{udot}n - (1/2*beta -1){u dotdot}n
		
    	delete[] VEL10;
    	delete[] ACC10;
    	delete[] DIS10;
    	delete[] DIS_plus10;

    	// Equation 11
    	double* ACC11 = new double[Nx*n];  		Matrix_Copy(ACC11, ACC, vec_size);
    	double* ACC_plus11 = new double[Nx*n];  Matrix_Copy(ACC_plus11, ACC_plus, vec_size);

    	vector_by_scalar(ACC11, delta_t*(1-gamma), vec_size);
    	vector_by_scalar(ACC_plus11, delta_t*gamma, vec_size);
    	adding_three_arrays(ACC_plus, VEL, ACC11, ACC_plus11, vec_size);

    	delete[] ACC11;
    	delete[] ACC_plus11;

    	if(i < 15){
    		std::cout << "\n\n\n ---------------------- Time step: " << i << " ---------------------" << std::endl;
    		std::cout << "Displacement: " << std::endl;
    		for(int i = 1; i < Nx*n;i+=3)
    			std::cout << DIS_plus[i] << "    ";
    	}

    	// ---------------------- Assigning for next loop --------------

    	Matrix_Copy(DIS, DIS_plus, vec_size);
    	Matrix_Copy(VEL, VEL_plus, vec_size);
    	Matrix_Copy(ACC, ACC_plus, vec_size);
	}


	std::cout << "Final Displacement: " << std::endl;
    for(int i = 1; i < Nx*n;i+=3)
    	std::cout << DIS_plus[i] << "    ";

	delete[] Fi;
	delete[] DIS;
	delete[] VEL;
	delete[] ACC;
	delete[] DIS_plus;
	delete[] VEL_plus;
	delete[] ACC_plus;

	return;
}