# include <iostream>
# include <iomanip>
# include <array>
# include <algorithm>
# include "math.h"
# include <string>

//# include "task1.h"
# include "task2.h"


int main(){

	// ----------------------- Task 1 (a) ----------------------------
	int Nx = 2;
    double b = 0.1;
    double h = 0.12;
    double L = 10.0;
    double A = 0.1*0.12;
    double I = b*pow(h, 3)/12;
    double E = 210000000000;
    double l = L/Nx;


    Nx = Nx + 1;
    // ----------------------- Variables needed ----------------------
    const double Qx = 0;       // no axial force in Task 1;    
    const double Qy = 10^3;    //  originally in N/mm   

    // ---------------------- Further Task 1 -------------------------
    T1_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy);

    // ---------------------- Further Task 2 -------------------------
    // need to add choice of static or dynamic
    double T = 1;
    double Nt =  5;
    double rho = 7850;

    T2_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, T, Nt, rho);

    std::cout << "\n\nEnd of program" << std::endl;
    return 0;
}


