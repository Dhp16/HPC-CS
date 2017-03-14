// bsv.cpp

// need to write function to write switched matrix

// need to add the force 

#include <fstream>


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
}

// ------------------------- Tools ----------------------------------
void output_array(const double *A, const std::string s, const int size){
    std::cout << "Printing: " << s << std::endl;
    for(int i = 0; i < size*size; i++){
        if(i != 0 && i % size == 0)
            std::cout << "\n";

        std::cout << std::setw(15);
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
        arr[i] = 0;
    return;
}

//  ------------------------ Task 1 (a) -----------------------------
double Value_Input(std::string name){
    std::string input;
    double value;
    do {
        std::cout << "Please enter desired "<< name <<": ";
        getline (std::cin, input);
        value = std::stod(input);
        if (value <= 0.0) {
            std::cout << name << " must be positive!" << std::endl;
        }
     } while (value <= 0.0);
    return value;
}
int Nx_Input(){
    std::string input;
    int value = 0;
    do {
        std::cout << "Please enter desired Nx: ";
        getline (std::cin, input);
        value = stoi(input);
        if (value <= 0.0) {
            std::cout << "Nx must be positive!" << std::endl;
        }
        else if(value % 2 == 0.){
            std::cout << "Nx must be impair." << std::endl;
        }
     } while (value <= 0.0 || value % 2 == 0.);
    return value;
}
// -------------------------- Task 1 (b) ---------------------------------
void pp_Fe(double Fe[6], double l, int Nx, double Qy, double Qx){      
    Fe[0] = Qx/2;
    Fe[1] = Qy/2;
    Fe[2] = Qy*l/12;
    Fe[3] = Qx/2;
    Fe[4] = Qy/2;
    Fe[5] = -Qy*l/12;
    for(int i = 0; i < 6; i++)                                     // could create an operator for this
        Fe[i] = Fe[i]*l;
    
    return; 
}
void pp_Ke(double* Ke, const double A,const double E,const double I,const double l){  // Populating the element stiffness matrix
    //double Ke[36] = {0};
    int n = 6;
    Ke[0] = A*E/l;
    Ke[1+n] = 12*E*I/pow(l, 3);
    Ke[2+2*n] = 4*E*I/l;
    Ke[3+3*n] = A*E/l;
    Ke[4+4*n] = Ke[1+n];
    Ke[5+5*n] = Ke[2+2*n];
    //
    Ke[0+3*n] = -A*E/l;
    Ke[1+2*n] = 6*E*I/pow(l,2);
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -6*E*I/pow(l,2);
    Ke[2+5*n] = 2*E*I/l;
    Ke[4+5*n] = -6*E*I/pow(l, 2);
    for(int i = 0; i < 6; i++){                         // Filling reamining fields knowing the matrix is symetric
        for(int j = 0; j < 6; j++){
            if(i == j)
                continue;
            if(i > j)
                Ke[i+j*n] = Ke[j+i*n];
        }
    }
    //output_array(Ke, "Ke", 6);
    return;
}
// ------------------------- Banded Non Symetric ------------------------------
void Global_Force_Vector(double* G_Fe, const double Fe[6], const int Nx, const int k, const double Fy){
    int n = 3;
    for(int i = 0; i < n*Nx; i = i+3){
        G_Fe[i] = Fe[0]+Fe[3];
        G_Fe[i+1] = Fe[1]+Fe[4];
        G_Fe[i+2] = Fe[2]+Fe[5];
    }
    int index = int(Nx*n/2);
    G_Fe[index] = G_Fe[index] + Fy;
    return;
}
void bfl(double* bfl, const double* Ke, const int Nx){                         // banded first line
    int nK = 6;
    int n = 3;
    int k = 4;
    for(int i = 0; i < Nx*3; i++){
        if(i < k){
            bfl[i] = 0;
            continue;
        }
        if((i+1) % 3 == 0)
            bfl[i] = Ke[5+nK];
        else
            bfl[i] = 0;
    }
}
void bsl(double* bsl, const double* Ke, const int Nx){
    int nK = 6;
    int n = 3;
    int k = 3;
    for(int i = 0; i < k; i++)
        bsl[i] = 0;
    for(int i = k; i < Nx*3; i = i+3){
        bsl[i] = Ke[3*nK];
        bsl[i+1] = Ke[1+4*nK];
        bsl[i+2] = Ke[2+5*nK];
    }
    return;
}
void btl(double* btl, const double* Ke, const int Nx){
    int nK = 6;
    int n = 3;
    int k = 2;
    for(int i = 0; i < k; i++)
        btl[i] = 0;
    for(int i = k; i < Nx*3; i++){
        if((i-1) % 3 == 0)
            btl[i] = Ke[2+4*nK];
        else 
            btl[i] = 0;
    }
}
void bfol(double* bfl,const double* Ke, const int Nx){
    for(int i = 0; i < Nx*3; i++)
        bfl[i] = 0;
    return;
}
void bmid(double* bmid, const double* Ke, const int Nx){
    int nK = 6;
    for(int i = 0; i < Nx*3; i = i+3){
        bmid[i] = 2*Ke[0];
        bmid[i+1] = 2*Ke[1+nK];
        bmid[i+2] = 2*Ke[2+2*nK]; 
    }
    return;
}
void Global_Stiffness_Matrix(double* G_Ke, const double* Ke, const int Nx, const int k){
    // k is passed through as the number of super diagonals, hence the number of 0 rows needed in the global stiffness matrix

    const int n = 3;
    double* fl = new double[Nx*n];
    double* se = new double[Nx*n];
    double* th = new double[Nx*n];
    double* fo = new double[Nx*n];
    double* mid = new double[Nx*n];

    bfl(fl, Ke, Nx);
    bsl(se, Ke, Nx);
    btl(th, Ke, Nx);
    bfol(fo, Ke, Nx);
    bmid(mid, Ke, Nx);

    int counter = 1;
    for(int i = 0; i < 9+k; i++){
        if(i == 0+k || i == 8+k){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = fl[j];
            }
        }
        else if(i == 1+k || i == 7+k){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = se[j];
            }
        }
        else if(i == 2+k || i == 6+k){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = th[j];
            }
        }
        else if(i == 3+k || i == 5+k){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = fo[j];
            }
        }
        else if(i == 4+k){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = mid[j];
            }
            back_array(fo, 1, Nx);
            back_array(th, 2, Nx);
            back_array(se, 3, Nx);
            back_array(fl, 4, Nx);
        }
        else if(i < k){                     // adding rows for multipliers
             for(int j = 0; j < Nx*n; j++)
                G_Ke[i*Nx*n+j] = 0;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////

// ----------------------- Task 1 (d) ------------------------------------
void write_to_file(double* Sol, const int n, const double L, const int Nx){
    std::ofstream myfile;
    std::string file_name = "T1_csol.txt";
    myfile.open(file_name);
    myfile << 0 << std::endl;
    for(int i = 1; i < n; i+=3){
        myfile << Sol[i] << std::endl;
    }
    myfile << 0 << std::endl;
    myfile.close();
    std::cout << "Solution written to file: " << file_name << std::endl;

    std::vector<double> x_locs;
    double delta_x = L/(Nx+1);
    for(int i = 0; i < Nx+2; i++){
        x_locs.push_back((double)i*delta_x);
    }
    std::ofstream myfile2;
    std::string file_name2 = "X_csol.txt";
    myfile2.open(file_name2);
    for(int i = 0; i < x_locs.size(); i++)
        myfile2 << x_locs[i] << std::endl;
    myfile2.close(); 

    return;
}
void Solve_Eq2(double* G_Fe, double* N, const int Nx, const double L){                                       // [K]{u} = {F}
    // DGBMV
    // http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_ga0dc187c15a47772440defe879d034888.html#ga0dc187c15a47772440defe879d034888
    // DGBSV
    // http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7
    int* ipiv = new int[Nx*3];

    int n = 3;
    int J = Nx*n;
    int kl = 4;
    int ku = 4;
    int lda = 1+2*kl+ku;
    int ldb = J;
    int nrhs = 1;
    int info = 16;

    F77NAME(dgbsv) (J, kl, ku, nrhs, N, lda, ipiv, G_Fe, ldb, &info);
    //F77NAME(dgbsv) (Nx*3, 4, 4, 1, N, Nx*3, ipiv, G_Fe, Nx*3, &info); 
    if(info == 0){
        std::cout << "Computation of the solution successful" << std::endl;
        write_to_file(G_Fe, J, L, Nx);
    }
    else
        std::cout << "Unable to solve system" << std::endl;

    delete[] ipiv;
    return;
}   
void Matrix_Transformer(double* N, const double* M, const int Nx, const int k){
    const int n = 3;

    for(int i = 0; i < Nx*n; i++){
        for(int j = 0; j < 9+k; j++){
            N[i*(9+k)+j] = M[j*Nx*n+i]; 
        }
    }

    /*std::cout << "Global stiffness Switched" << std::endl;
    for(int i = 0; i < Nx*n*13; i++){
        if(i % 13 == 0)
            std::cout << std::endl;
        std::cout << std::setw(5);
        std::cout << N[i] << " " ;
    }*/

    return;
}
// ----------------------- Task 1 (e) ------------------------------------
void Analytical_Sol(const int Nx, const double L, const double Qy, const double E, const double I, const double Fy){
    std::vector<double> x_locs;
    std::vector<double> asol;

    double delta_x = L/(100);
    for(int i = 0; i < 100; i++){
        x_locs.push_back((double)i*delta_x);
    }

    for(int i = 0; i < x_locs.size()/2+1; i++){
        double x = x_locs[i];
        double u_distributed = Qy*pow(x, 2)*pow((L-x), 2)/(24*E*I);
        double u_concentrated = Fy*pow(x, 2)*(3*L-4*x)/(48*E*I);
        double u_total = u_distributed + u_concentrated;
        asol.push_back(u_total);
    }
    int k = asol.size();
    for(unsigned int i = 0; i < x_locs.size()-k; i++){
        asol.push_back(asol[k-2-i]);
    }

    std::ofstream myfile;
    std::string file_name = "T1_asol.txt";
    myfile.open(file_name);
    for(int i = 0; i < asol.size(); i++)
        myfile << asol[i] << std::endl;
    myfile.close(); 

    std::ofstream myfile2;
    std::string file_name2 = "X_asol.txt";
    myfile2.open("X.txt");
    for(int i = 0; i < x_locs.size(); i++)
        myfile2 << x_locs[i] << std::endl;
    myfile2.close(); 

    //std::cout << "Analytical Solution written to: " << file_name << std::endl;
    return;
}


// -----------------------------------------------------------------------
void T1_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy){
    //int Nx = Nx_Input();
    //double A = Value_Input("A");
    //double I = Value_Input("I");
    //double E = Value_Input("E");
    //double L = Value_Input("L");
    // ----------------------- Task 1 (a) ----------------------------
    const int n = 3;
    const int k = 4;
    // ----------------------- T1b -----------------------------------
    double Ke[36] = {0};
    double Fe[6] = {0};
    pp_Fe(Fe, l, Nx, Qy, Qx);
    pp_Ke(Ke, A, E, I, l);
    // ----------------------- Task 1 (c) ----------------------------
    //double* G_Fe = new double[Nx*6];
    //double* G_Ke = new double[Nx*6*Nx*6];
    //Global_Fe(G_Fe, l, Nx, Qy);
    //Global_Ke(Ke, Nx, A, E, I, l);

    //double* G_Fe_nBC = new double[(Nx-1)*6];                // without the boundary conditions
    //double* K_Fe_nBC = new double[(Nx-1)*6*(Nx-1)*6];       // without the boundary conditions
    //GFe_Take_Off_BCs(G_Fe, G_Fe_nBC, Nx);

    double* G_Ke = new double[Nx*n*(10+k)];
    double* N = new double[Nx*n*(10+k)];
    double* G_Fe = new double[n*Nx+k];

    Global_Force_Vector(G_Fe, Fe, Nx, k, Fy);
    Global_Stiffness_Matrix(G_Ke, Ke, Nx, k);

    Matrix_Transformer(N, G_Ke, Nx, k);
    Solve_Eq2(G_Fe, N, Nx, L);
    Analytical_Sol(Nx, L, Qy, E, I, Fy);

    std::cout << "\nEnding task 1" << std::endl;

    delete[] N;
    delete[] G_Ke;
    delete[] G_Fe;

    return;
}

// Nx number of elements, number of nodes  = number of elements + 1