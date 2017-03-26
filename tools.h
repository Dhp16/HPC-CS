// This file is used to store functions that are used for multiple tasks



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
void Inverse_Diagonal_Matrix(double* M, const int Nx){
    // M any diagonal matrix
    // inverse of a diagonal matrix is obtained by replacing each element by its inverse
    int n = 3;
    for(int i = 0; i < Nx*n; i++){
        M[i] = 1./M[i];
    }
    return;
}
void Diagonal_by_vector_Imp(double* diagonal, const double* vector, const int Nx){
    const int n = 3;
    for(unsigned int i = 0; i < Nx*n; i++){
        diagonal[i] = diagonal[i]*vector[i];
    }
    return;
}
void Matrix_by_Scalar(double* M, const double u, const int Nx){
    const int n = 3;
    for(int i = 0; i < Nx*n; i++){
        M[i] = M[i]*u;
    }
    return;
}
void Matrix_Add_Diagonal(double* AA, const double* M, const int Nx, const int k = 4){
    // k: number of upper diagonals
    int n = 3;
    int counter = 0;
    for(int i = k*Nx*n; i < (k+1)*Nx*n; i++){
        AA[i] = AA[i] + M[counter]; 
        counter++;
    }
    return;
}
void Matrix_Copy(double* copy, const double* original, const int size){
    for(int i = 0; i < size; i++)
        copy[i] = original[i];
}
void print(double* A, const int size){
    for(int i = 0; i < size; i++)
        std::cout << std::setw(12) << A[i] << "   ";
    std::cout << std::endl;
}
void output_array(const double *A, const std::string s, const int cols, const int lines){
    std::cout << "Printing: " << s << std::endl;
    for(int i = 0; i < cols*lines; i++){
        if(i != 0 && i % cols == 0)
            std::cout << "\n";
        std::cout << std::setw(12);
        std::cout << A[i] << " ";
    }
    std::cout << "\n";
}
void back_array(double* arr, const int k, const int Nx){
    const int n = 3;
    double* copy = new double[Nx*n];
    copy = arr;
    for(int i = k; i < Nx*n; i++){
        arr[i-k] = copy[i];
    }
    for(int i = Nx*n-k; i <Nx*n; i++){
        arr[i] = 0;
    }
    return;
}
void initialise_dynamic_array(double* arr, const int Nx, const int n, const int m){
    for(int i = 0; i < pow(Nx*n, m); i++)
        arr[i] = 0.0;
    return;
}
void initialise_dynamic_array(double* arr, const int size){
    for(int i = 0; i < size; i++)
        arr[i] = 0.0;
    return;
}
void Matrix_Transformer(double* N, const double* M, const int Nx, const int k = 0){
    const int n = 3;

    for(int i = 0; i < Nx*n; i++){
        for(int j = 0; j < 9+k; j++){
            N[i*(9+k)+j] = M[j*Nx*n+i]; 
        }
    }

    return;
}
void Matrix_Transformer_Imp(double* M, const int Nx, const int k = 0){              // Number of rows of 0s
    const int n = 3;
    double* X = new double[Nx*n*(9+k)];
    for(int i = 0; i < Nx*n*(9+k); i++)
        X[i] = M[i];

    for(int i = 0; i < Nx*n; i++){
        for(int j = 0; j < 9+k; j++){
            M[i*(9+k)+j] = X[j*Nx*n+i]; 
        }
    }
    return;
}
void vector_by_scalar(double* vector, const double scalar, const int size){
    for(int i = 0; i < size; i++)
        vector[i] = vector[i]*scalar;
    return;
}
void adding_three_arrays(double* result, const double* A, const double* B, const double* C, const int size){
    // /!\ used only for TASK 2
    for(int i = 0; i < size; i++){
        result[i] = A[i] - B[i] + C[i];
    }
}
void adding_to_array(double* A, const double* B, const int size){
    for(int i = 0; i < size; i++){
        A[i] = A[i] + B[i];
    }
}
void sub_from_vector(double* A, const double* B, const int size){
    for(int i = 0; i < size; i++){
        A[i] = A[i] - B[i];
    }
}
void Diagonal_by_Vector(double* result, const double* diagonal, const double* vector, const int size){
    for(unsigned int i = 0; i < size; i++){
        result[i] = diagonal[i]*vector[i];
    }
    return;
}
void write_task_solution(const double* Sol, const int size, const std::string task){
    std::ofstream file;
    std::string file_name = "Task" + task + "_Solution.txt"; 
    file.open(file_name);
    for(int i = 1; i < size; i+=3)
        file << Sol[i] << std::endl; 
    file.close();
    std::cout << "Vertical displacement of the beam written to: " << file_name << std::endl;
    return;
}