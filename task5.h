/*  
 Imperial College London
 HPC Assignment Task 5
 Dominic Pickford 
 01272723

 This is the source file declaring all the functions that are key to task 4

 The main loop is very close to that designed for task 2, with the addition of communication between the processes

 Inspiration from example pbsv.cpp provided
*/

#ifndef TASK5_H_INCLUDED
#define TASK5_H_INCLUDED

void T5_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho);

#endif 