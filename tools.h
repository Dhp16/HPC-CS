// This file is used to store functions that are used for multiple tasks

#ifndef TOOLS_H_INCLUDED
#define TOOLS_H_INCLUDED

# include <vector>
# include <string>
# include "math.h"
# include <iostream>
# include <fstream>

#define F77NAME(x) x##_
extern "C" {
    void F77NAME(dgbsv)(const int& n, 
                    const int& kl, 
                    const int& ku, 
                    const int& nrhs, 
                    const double * A, 
                    const int& ldab, 
                    int * ipiv, double * B, 
                    const int& ldb, int* info);

    void F77NAME(dgbmv)(const char& trans, 
                    const int& m, 
                    const int& n,
                    const int& kl, 
                    const int& ku,
                    const double& alpha, const double* a, 
                    const int& lda,
                    const double* x, const int& incx, 
                    const double& beta,
                    double* y, const int& incy);

    void F77NAME(dpbsv)(const char &UPLO,
                    const int& N,
                    const int& KD,
                    const int& NRHS,
                    double* A,
                    const int& LDA,
                    double* B,
                    const int LDB,
                    int& INFO);
    void F77NAME(pdgbsv)(const int& n, const int& kl, const int& ku, 
                    const int& nrhs, const double * A, const int& ja,
                    const int* desca, int * ipiv, double * B, const int& ib,
                    const int* descb, double* work, const int& lwork, 
                    int* info); 
    void Cblacs_get(int, int, int*);
    void Cblacs_pinfo(int*, int*);
    void Cblacs_gridinit(int*, char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_exit(int);
    void Cblacs_gridexit(int);

    //void F77NAME(dpbsv)(const char& uplo,const int& n,const int& ku,const int& nrhs,double *A,const int& ldab,double* B,const int& ldb, int& info);
}  

// ------------------------- Tools ----------------------------------
void Inverse_Diagonal_Matrix(double* M, const int Nx);
// Multiply a diagonal matrix by a vector: element by element
void Diagonal_by_vector_Imp(double* diagonal, const double* vector, const int Nx);
// Multiply every element in a matrix by a scalar
void Matrix_by_Scalar(double* M, const double u, const int Nx);
// Add a vector to the diagonal of a matrix (already banded)
void Matrix_Add_Diagonal(double* AA, const double* M, const int Nx, const int k);
// Copy a Matrix 
void Matrix_Copy(double* copy, const double* original, const int size);
// alternative function to print a matrix to the console (capable of handling more matrices)
void output_array(const double *A, const std::string s, const int cols, const int lines);
// shift the values of an array forward by a given number of columns
void back_array(double* arr, const int k, const int Nx);
// set all 0s to a dynamic array
void initialise_dynamic_array(double* arr, const int Nx, const int n, const int m);
// set all 0s to a dynamic array
void initialise_dynamic_array(double* arr, const int size);
// Shift a matrix from Row Major to Column Major, N output array (CM), M input, Nx number of elements, k number of 0 lines (included for use of BLAS/LAPACK functions) 
void Matrix_Transformer(double* N, const double* M, const int Nx, const int k);
// Same as above but replaces the array rather than requiring another
void Matrix_Transformer_Imp(double* M, const int Nx, const int k);
// Multiply a vector by a scalar, element by element
void vector_by_scalar(double* vector, const double scalar, const int size);
// Function used specifically for Task 2, combining three arrays
void adding_three_arrays(double* result, const double* A, const double* B, const double* C, const int size);
// adding two arrays together element by element
void adding_to_array(double* A, const double* B, const int size);
// substract a constant array 
void sub_from_vector(double* A, const double* B, const int size);
void Diagonal_by_Vector(double* result, const double* diagonal, const double* vector, const int size);
void write_task_solution(const double* Sol, const int size, const std::string task);

#endif