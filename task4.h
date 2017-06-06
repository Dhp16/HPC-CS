/* 
 Imperial College London
 HPC Assignment Task 4
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 4

 The main loop is very close to that designed for task 2, with the addition of communication between the processes
*/


// Parallelise your program

/* each process has unique ID --> rank
 each works on different part of the problem


 link with MPI libraries or compile with mpicxx

 run with mpiexec command

	mpicxx hello.cpp -o hello
	mpiexcec -np 2 ./hello

each message contains:
- sender rank
- destination rank
- message datatype
- message size
- message location
- message tag

 */ 

//void PARA_build_M(const double rho, const double A, const double l, const int PARA_size, const int rank){

// /} 


#ifndef TASK4_H_INCLUDED
#define TASK4_H_INCLUDED

void write_task4(const int Nx, double* Solution, const int vec_size, const double l);

void PARA_Force_Vector(double* PARA_G_Fe, const double l, const int PARA_Nx, const int Nx, const double Fy, const double Qy, const double Qx, const int rank, const double scalar);

void T4_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
	const double T, const double Nt, const double rho);


#endif 