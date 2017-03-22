/* 
 Imperial College London
 HPC Assignment Task 3
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 3
*/


void print(double* A, const int size){
	for(int i = 0; i < size; i++)
		std::cout << std::setw(12) << A[i] << "   ";
	std::cout << std::endl;
}

void K_Effective(double* K_effective, const double* G_Mm, const int Nx, const double delta_t){
	// assuming the K_effective that is input is already equal to G_Ke
	int n = 3;
	double* G_Mm_copy = new double[Nx*n];
	Matrix_Copy(G_Mm_copy, G_Mm, Nx*n);  

	double u = 1/(0.25*pow(delta_t,2));

	Matrix_by_Scalar(G_Mm_copy, u, Nx);
	Matrix_Add_Diagonal(K_effective, G_Mm_copy, Nx, 8);
	delete[] G_Mm_copy;
}

void write_task3(std::string file_name, std::vector<double> tsteps, std::vector<double> dispMid){
	std::ofstream file;
	file_name = "Task3_" + file_name + ".txt"; 
	file.open(file_name);
	for(int i = 0; i < dispMid.size(); i++)
		file << tsteps[i] << std::endl; 
	for(int i = 0; i < dispMid.size(); i++)
		file << dispMid[i] << std::endl; 
	file.close();
	return;
}



void T3_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho){

	const double beta = 0.250; 
	const double gamma = 0.5;

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

	double* K_effective = new double[vec_size*13];
	Global_Stiffness_Matrix(K_effective, Ke, Nx, 4);

	K_Effective(K_effective, G_Mm, Nx, delta_t);		// inside the loop beacuse DGBSV mangles with it

	Matrix_Transformer_Imp(K_effective, Nx, 4);	

	std::vector<double> tsteps;
	std::vector<double> dispMid;

	double* DIS = new double[vec_size];
	double* VEL = new double[vec_size];
	double* ACC = new double[vec_size];
	initialise_dynamic_array(DIS, Nx, n, 1);
	initialise_dynamic_array(VEL, Nx, n, 1);
	initialise_dynamic_array(ACC, Nx, n, 1);
	double* DIS_plus = new double[vec_size];
	double* VEL_plus = new double[vec_size];
	double* ACC_plus = new double[vec_size];
	initialise_dynamic_array(DIS_plus, Nx, n, 1);
	initialise_dynamic_array(VEL_plus, Nx, n, 1);
	initialise_dynamic_array(ACC_plus, Nx, n, 1);


	for(int i = 0; i < Nt; i++){
		if(i < 5){std::cout << "\n\n\n ---------------------- Time step: " << i << " ---------------------" << std::endl;}
		double t_now = i*delta_t;

		double* Fi = new double[vec_size];
		initialise_dynamic_array(Fi, Nx, n, 1);
		if(2*(t_now+delta_t)/T > 1)
			Global_Force_Vector(Fi, Fe, Nx, k, Fy, 1);
		else
			Global_Force_Vector(Fi, Fe, Nx, k, Fy, 2*(t_now+delta_t)/T);

		if(i < 5){std::cout << "Force" << std::endl; print(Fi, vec_size);}

		// -------------------- Equation 12	---------------------
		double* A = new double[vec_size];
		for(int i = 0; i < vec_size; i++)
			A[i] = (1/(beta*pow(delta_t, 2))) * DIS[i] + (1/(beta*delta_t)) * VEL[i] + (0.5/beta-1)*ACC[i];		

		if(i < 5){std::cout << "A" << std::endl; print(A, vec_size);}

		double* B = new double[vec_size];
		for(int i = 0; i < vec_size; i++)
			B[i] = Fi[i] + A[i]*G_Mm[i];																	// B = {F}n+1  + [M]{......}
		
		if(i < 5){
		std::cout << "Global Mass Matrix" << std::endl; print(G_Mm, vec_size);			
		std::cout << "B" << std::endl; print(B, vec_size);}

		delete[] A;												

		int* ipiv = new int[Nx*3];
    	int info = 1;

    	double* K_effective_copy = new double[vec_size*13];
    	Matrix_Copy(K_effective_copy, K_effective, vec_size*13); // protecting K_effective from the changes incureed by dgbsv
    	
     	F77NAME(dgbsv) (vec_size, 4, 4, 1, K_effective_copy, 1+2*4+4, ipiv, B, vec_size, &info);	// result should be ported to DIS_plus10--> save KEff for next loop


    	if(i == 0 && info != 0)
    		std::cout << "dgbsv failed!" << std::endl;

    	Matrix_Copy(DIS_plus, B, vec_size);											// succesfully copies	
    	
    	if(i < 5){std::cout << "DIS_Plus" << std::endl; print(B, vec_size);}
    	delete[] K_effective_copy;
    	delete[] B;


		// ---------------------- Equation 10 ------------------- 
    	for(int i = 0; i < vec_size; i++)
    		ACC_plus[i] = ((1/(beta*pow(delta_t, 2)))*(DIS_plus[i] - DIS[i])) - (1/(beta*delta_t))*VEL[i] - ((1/(2*beta)-1)*ACC[i]);   

		if(i < 5){std::cout << "ACC_plus" << std::endl; print(ACC_plus, vec_size);}
    	// ----------------------- Equation 11 -------------------
    	for(int i = 0; i < vec_size; i++)
    		VEL_plus[i] =       VEL[i]      +         (delta_t*(1-gamma))*ACC[i]     +    delta_t*gamma*ACC_plus[i];
    	if(i < 5){std::cout << "VEL_plus" << std::endl; print(VEL_plus, vec_size);}

    	if(i % 1000 == 0){
    		std::cout << "\n\n\n ---------------------- Time step: " << i << " ---------------------" << std::endl; 
    		std::cout << "Force" << std::endl;
    		for(int i = 0; i < vec_size; i++)
    			std::cout << Fi[i] << "   ";
    		std::cout << "\nDisplacement: " << std::endl;
    		for(int i = 1; i < Nx*n;i+=3)
    			std::cout << DIS_plus[i] << "    ";
    		std::cout << "\nVelocity: " << std::endl;
    		for(int i = 1; i < Nx*n;i+=3)
    			std::cout << VEL_plus[i] << "    ";
    		std::cout << "\nAcceleration: " << std::endl;
    		for(int i = 1; i < Nx*n;i+=3)
    			std::cout << ACC_plus[i] << "    "; 
    	}

    	// ---------------------- PUSH BACK FOR PLOTS ------------------
    	tsteps.push_back(t_now);
    	dispMid.push_back(DIS_plus[(Nx*3/2)-2]);

    	// ---------------------- Clearing Memory ----------------------
    	delete[] Fi;
    	// ---------------------- Assigning for next loop --------------
    	Matrix_Copy(DIS, DIS_plus, vec_size);
    	Matrix_Copy(VEL, VEL_plus, vec_size);
    	Matrix_Copy(ACC, ACC_plus, vec_size);
    	initialise_dynamic_array(DIS_plus, Nx, n, 1);	// making sure these are returned to 0 before next iteration
		initialise_dynamic_array(VEL_plus, Nx, n, 1);
		initialise_dynamic_array(ACC_plus, Nx, n, 1);

	}

	std::string file_name = "HT";
	write_task3(file_name, tsteps, dispMid);


	std::cout << "\n\n\n----- Final Displacement: -----" << std::endl;
    for(int i = 1; i < Nx*n;i+=3)
    	std::cout << DIS_plus[i] << "    ";

	delete[] DIS;
	delete[] VEL;
	delete[] ACC;
	delete[] DIS_plus;
	delete[] VEL_plus;
	delete[] ACC_plus;

	return;
}