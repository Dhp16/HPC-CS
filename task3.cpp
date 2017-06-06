/* 
 Imperial College London
 HPC Assignment Task 3
 Dominic Pickford 
 01272723

 This is the source file declaring all the functions that are key to task 3

 DGBSV documentation
 http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7
*/

# include "task3.h"
# include "tools.h"
# include "task1.h"
# include "task2.h"
# include "math.h"

void K_Effective(double* K_effective, const double* G_Mm, const int Nx, const double delta_t){
	// assuming the K_effective that is input is already equal to G_Ke
	int n = 3;											// columns per node
	double* G_Mm_copy = new double[Nx*n];
	Matrix_Copy(G_Mm_copy, G_Mm, Nx*n);  

	double u = 1/(0.25*pow(delta_t,2));					// define the time factor

	Matrix_by_Scalar(G_Mm_copy, u, Nx);					// pre-multiply G_Mm by the time factor
	Matrix_Add_Diagonal(K_effective, G_Mm_copy, Nx, 8);	// Add G_Mm diagonal to K_effective
	delete[] G_Mm_copy;									// clear memory
}
void write_task3_displacements(std::vector<double> tsteps, std::vector<double> dispMid){	
	std::ofstream file;
	std::string file_name = "Task3_Mid_point_vs_Time.txt"; 
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

	std::cout << "\n---------------------- Task 3 ------------------------" << std::endl;

	const double beta = 0.250; 									
	const double gamma = 0.5;

	int n = 3;
	int vec_size = Nx*n; 
	double delta_t = T/Nt;										// time step 

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

	//	DIS short for displacement     DIS_plus = {u}n+1
	//	VEL short for velocity         VEL_plus = {udot}n+1
	//	ACC short for acceleration     ACC_plus = {udotdot}n+1


	double* DIS = new double[vec_size];				
	double* VEL = new double[vec_size];
	double* ACC = new double[vec_size];
	initialise_dynamic_array(DIS, Nx, n, 1);	// Set to 0
	initialise_dynamic_array(VEL, Nx, n, 1);
	initialise_dynamic_array(ACC, Nx, n, 1);
	double* DIS_plus = new double[vec_size];
	double* VEL_plus = new double[vec_size];
	double* ACC_plus = new double[vec_size];
	initialise_dynamic_array(DIS_plus, Nx, n, 1);
	initialise_dynamic_array(VEL_plus, Nx, n, 1);
	initialise_dynamic_array(ACC_plus, Nx, n, 1);

	std::vector<double> tsteps;
	std::vector<double> dispMid;  

	// Control of loading time
	/* 
		As it stands the load will be applied linearly from [0; T], however, changing the "factor" below to a higher number will decrease the loading time. A vector of 4000 different 
		factors were used to plot the results shown for amplitude against loading time in the report.
	*/
	double factor = 1.0;

	bool loading_finished = false;
	double loading_time = 0;
	double end_loading_time_step = 0;
	

	for(int i = 0; i < Nt; i++){
		double t_now = i*delta_t;																			// real time in this timestep

		double* Fi = new double[vec_size];
		initialise_dynamic_array(Fi, vec_size);																// initilise force vector to 0
		if(factor*(t_now+delta_t)/T >= 1)																	// if the maximum loading has been reached, keep scalling cst at 1
			Global_Force_Vector(Fi, Fe, Nx, k, Fy, 1);
		else
			Global_Force_Vector(Fi, Fe, Nx, k, Fy, factor*(t_now+delta_t)/T);								// if the maximum loading is not yet reached scale by: 2*(t_now+delta_t)/T

		// -------------------------------------------------------
		 if(factor*t_now/T >= 1.0 && loading_finished == false){											// flag the fact the maximum loading has been reached
 			loading_finished = true;
 			loading_time = t_now; 
 			end_loading_time_step = i; 
 		} 
		// -------------------------------------------------------

		// -------------------- Equation 12	---------------------
		double* A = new double[vec_size];																	
		for(int i = 0; i < vec_size; i++)
			A[i] = (1/(beta*pow(delta_t, 2))) * DIS[i] + (1/(beta*delta_t)) * VEL[i] + (0.5/beta-1)*ACC[i];		

		double* B = new double[vec_size];
		for(int i = 0; i < vec_size; i++)
			B[i] = Fi[i] + A[i]*G_Mm[i];																	// B = {F}n+1  + [M]{......}

		delete[] A;												

		int* ipiv = new int[Nx*3];
    	int  info = 1;																				// default info to 1 in case there is a failure

    	double* K_effective_copy = new double[vec_size*13];
    	Matrix_Copy(K_effective_copy, K_effective, vec_size*13); 									// protecting K_effective from the changes incureed by dgbsv
    	
     	F77NAME(dgbsv) (vec_size, 4, 4, 1, K_effective_copy, 1+2*4+4, ipiv, B, vec_size, &info);	// result should be ported to DIS_plus10--> save KEff for next loop
     	if(info != 0)
     		std::cout << "NOTICE: dgbmv info: " << info << std::endl; 

    	Matrix_Copy(DIS_plus, B, vec_size);															// copy B into Displacement[n+1] and delete B
    	
    	// Clear memory after DGBSV
    	delete[] K_effective_copy;
    	delete[] B;

		// ---------------------- Equation 10 ------------------- 
    	for(int i = 0; i < vec_size; i++)
    		ACC_plus[i] = ((1/(beta*pow(delta_t, 2)))*(DIS_plus[i] - DIS[i])) - (1/(beta*delta_t))*VEL[i] - ((1/(2*beta)-1)*ACC[i]);   

    	// ----------------------- Equation 11 -------------------
    	for(int i = 0; i < vec_size; i++)
    		VEL_plus[i] =       VEL[i]      +         (delta_t*(1-gamma))*ACC[i]     +    delta_t*gamma*ACC_plus[i];

    	// ---------------------- PUSH BACK FOR PLOTS ------------------
    	tsteps.push_back(t_now);
    	dispMid.push_back(DIS_plus[(Nx*3)/2]);

    	// ---------------------- Clearing Memory ----------------------
    	delete[] Fi;
    	// ---------------------- Assigning for next loop --------------
    	Matrix_Copy(DIS, DIS_plus, vec_size);
    	Matrix_Copy(VEL, VEL_plus, vec_size);
    	Matrix_Copy(ACC, ACC_plus, vec_size);
    	initialise_dynamic_array(DIS_plus, vec_size);	// making sure these are returned to 0 before next iteration
		initialise_dynamic_array(VEL_plus, vec_size);
		initialise_dynamic_array(ACC_plus, vec_size);
	}
	// ----------------------------------------------

    write_task_solution(DIS, vec_size, "3");			// write proof of solution solving (task 3 (a))
    //write_task3_displacements(tsteps, dispMid);			// wrire out displacement of mid points against time for verification

    // -------------------------- Clear Memory ------------------
    delete[] DIS;
	delete[] VEL;
	delete[] ACC;
	delete[] DIS_plus;
	delete[] VEL_plus;
	delete[] ACC_plus;


	std::cout << "------------------- End of Task 3 ----------------------" << std::endl;

	return;
}
