/* 
 Imperial College London
 HPC Assignment Task 2
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 2

 The variables used are all explained in function T2_inputs

 DGBMV documentation
 http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
*/

#ifndef TASK2_H_INCLUDED
#define TASK2_H_INCLUDED

# include <vector>

// ---------------------------------- Task 2 (b) ------------------
// Build Global Mass Matrix
void Global_Mass_Matrix(double* G_Mm, const int Nx, const double rho, const double A, const double l);
// -------------------------------- Task 2 (c) -------------------
// COEF1 = delta_t^2/[M]
void Build_Coef1(double* COEF1, const double* G_Mm, const int Nx, const double delta_t);
// COEF2 = [K] - 2/delta_t^2 *[M]
void Build_Coef2(double* COEF2, const double* G_Ke, const double* G_Mm, const int Nx, const double delta_t);
// COEF3 = -1/deltat^2 * [M]
void Build_Coef3(double* COEF3, const double* G_Mm, const int Nx, const double delta_t);
// ------------------------------- Task 2 (e) -------------------- 
// Writting mid point deflection @ each time step
void write_task2_displacements(std::vector<double> tsteps, std::vector<double> dispMid);
// ------------------------------- Task 2 (f) --------------------
// Writting amplitude of mid-point oscillation for each loading time
void write_oscillations(const std::vector<double> amplitudes, const std::vector<double> lt);

// -------------------------- Primary Function -------------------
void T2_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho);

/* 
	As it stands the load will be applied linearly from [0; T], however, changing the "factor" below to a higher number will decrease the loading time. A vector of 4000 different 
	factors were used to plot the results shown for amplitude against loading time in the report.
*/

#endif




