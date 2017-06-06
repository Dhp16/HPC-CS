/* 
 Imperial College London
 HPC Assignment Task 3
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 3

 DGBSV documentation
 http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7
*/


#ifndef TASK3_H_INCLUDED
#define TASK3_H_INCLUDED


# include <vector>

void K_Effective(double* K_effective, const double* G_Mm, const int Nx, const double delta_t);
void write_task3_displacements(std::vector<double> tsteps, std::vector<double> dispMid);

void T3_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho);

#endif


