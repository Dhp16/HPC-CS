/* 
 Imperial College London
 HPC Assignment Task 1
 Dominic Pickford 
 01272723

 This is the source file declaring all the functions that are key to task 1

 DGBSV Documentation
 http://www.netlib.org/lapack/explore-html/d3/d49/group__double_g_bsolve_gafa35ce1d7865b80563bbed6317050ad7.html#gafa35ce1d7865b80563bbed6317050ad7

 */

# include "task1.h"
# include "tools.h"
# include "math.h"
# include <iostream>
# include <vector>
# include <fstream>

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
    return;
}
// ------------------------- Task 1 (c) ------------------------------
void Global_Force_Vector(double* G_Fe, const double Fe[6], const int Nx, const int k, const double Fy, const double scalar){
    int n = 3;
    for(int i = 0; i < n*Nx; i = i+3){
        G_Fe[i] = Fe[0]+Fe[3];
        G_Fe[i+1] = (Fe[1]+Fe[4])*scalar;
        G_Fe[i+2] = (Fe[2]+Fe[5])*scalar;
    }
    int index = int(Nx*n/2);
    G_Fe[index] = G_Fe[index] + Fy*scalar;

    return;
}

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

void bfl(double* bfl, const double* Ke, const int Nx){
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
    // One the 4 super diagonals and the mid line have been passed, the function "back_array" is used to 
    // offset the indices in the first to fourth line so that these line up correctly in the banded storage


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
                                                
    return;
}

// ----------------------- Task 1 (d) ------------------------------------
void write_to_file(double* Sol, const int n, const double L, const int Nx){
    std::ofstream myfile;
    std::string file_name = "Task1_Solution.txt";
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

    int* ipiv = new int[Nx*3];

    int n = 3;                 
    int J = Nx*n;               
    int kl = 4;
    int ku = 4;
    int lda = 1+2*kl+ku;
    int ldb = Nx*n;
    int nrhs = 1;
    int info = 16;

    F77NAME(dgbsv) (J, kl, ku, nrhs, N, lda, ipiv, G_Fe, ldb, &info);                                       

    if(info == 0){
        write_to_file(G_Fe, J, L, Nx);
    }
    else
        std::cout << "NOTICE: dgbsv info: " << info << ". Unable to solve system" << std::endl;

    delete[] ipiv;
    return;
}    
// ----------------------- Task 1 (e) ------------------------------------
void Analytical_Sol(const int Nx, const double L, const double Qy, const double E, const double I, const double Fy){
    std::vector<double> x_locs;                                 // x indices
    std::vector<double> asol;                                   // analytical solution

    double delta_x = L/(100);
    for(int i = 0; i < 100; i++){
        x_locs.push_back((double)i*delta_x);                                
    }
    for(int i = 0; i < x_locs.size()/2+1; i++){                         
        double x = x_locs[i];
        double u_distributed = Qy*pow(x, 2)*pow(L-x, 2)/(24*E*I);
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

    return;
}
// -----------------------------------------------------------------------
void T1_inputs(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy){
    std::cout << "\n\n---------------------- Task 1 ------------------------" << std::endl;
    
    // ----------------------- Task 1 (a) ----------------------------
    const int n = 3;
    const int k = 4;
    
    // ----------------------- Task 1 (b) ----------------------------
    double Ke[36] = {0};                    // declare elemental stiffness matrix
    double Fe[6] = {0};                     // declare elemental force vector
    pp_Ke(Ke, A, E, I, l);                  // populate elemental stiffness matrix
    pp_Fe(Fe, l, Nx, Qy, Qx);               // populate elemental force vector
    
    // ----------------------- Task 1 (c) ----------------------------
    double* G_Fe = new double[n*Nx+k];              // declare global force vector
    Global_Force_Vector(G_Fe, Fe, Nx, k, Fy);       // populate global force vector

    double* G_Ke = new double[Nx*n*(10+k)];         // declare global stiffness matrix
    double* N = new double[Nx*n*(10+k)];        
    Global_Stiffness_Matrix(G_Ke, Ke, Nx, k);       // populate glbal stiffness matrix
    Matrix_Transformer(N, G_Ke, Nx, k);             // switch to column major
    Solve_Eq2(G_Fe, N, Nx, L);                      // solve equation
    delete[] G_Ke;
    delete[] N;

    
    Analytical_Sol(Nx, L, Qy, E, I, Fy);
    std::cout << "------------------- Ending Task 1 ---------------------" << std::endl;

    delete[] G_Fe;

    return;
}
