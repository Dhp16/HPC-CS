/* 
 Imperial College London
 HPC Assignment Task 1
 Dominic Pickford 
 01272723

 This is the header file declaring all the functions that are key to task 1

 DGBSV Documentation
 http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7

 */

#ifndef TASK1_H_INCLUDED
#define TASK1_H_INCLUDED

#include <string>

//  ------------------------ Task 1 (a) -----------------------------
double Value_Input(std::string name);
int Nx_Input();
// -------------------------- Task 1 (b) ---------------------------------
void pp_Fe(double Fe[6], double l, int Nx, double Qy, double Qx);
void pp_Ke(double* Ke, const double A,const double E,const double I,const double l);
// ------------------------- Task 1 (c) ------------------------------
void Global_Force_Vector(double* G_Fe, const double Fe[6], const int Nx, const int k, const double Fy, const double scalar = 1.0);

/*
    Each of the functions below creates a line for the Global Stiffness matrix when banded from the different components of the elemental stiffness matrix
    bfl: banded first line
    bsl: banded second line
    btl: banded third line
    bfol: banded forth line
    bmid: banded mid line (diagonal) 
    These are then assembled together in Global_Stiffness_Matrix

    In each of the these functions, k denotes the indentation and n the number of columns per node 
*/

void bfl(double* bfl, const double* Ke, const int Nx);
void bsl(double* bsl, const double* Ke, const int Nx);
void btl(double* btl, const double* Ke, const int Nx);
void bfol(double* bfl,const double* Ke, const int Nx);
void bmid(double* bmid, const double* Ke, const int Nx);
void Global_Stiffness_Matrix(double* G_Ke, const double* Ke, const int Nx, const int k = 0);

// ----------------------- Task 1 (d) ------------------------------------
void write_to_file(double* Sol, const int n, const double L, const int Nx);
void Solve_Eq2(double* G_Fe, double* N, const int Nx, const double L);
// ----------------------- Task 1 (e) ------------------------------------
void Analytical_Sol(const int Nx, const double L, const double Qy, const double E, const double I, const double Fy);
// -----------------------------------------------------------------------
void T1_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, 
    const double Qx, const double Qy, const double Fy);


#endif