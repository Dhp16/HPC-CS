# include <iostream>
# include <iomanip>
# include <array>
# include <algorithm>
# include "math.h"
# include <string>

//# include "task1.h"
# include "task2.h"


void Input_Parameters_to_file(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy){

}

int main(){

	// ----------------------- Task 1 (a) ----------------------------
	int Nx = 4;
    double b = 0.10;
    double h = 0.120;
    double L = 10.0;
    double A = b*h;
    double I = b*pow(h, 3)/12;
    double E = 210000000000;
    double l = L/Nx;

    // ----------------------- Variables needed ----------------------
    const double Qx = 0.0;       // no axial force in Task 1;    
    const double Qy = -1000.0;    //  originally in N/mm   
    const double Fy = -1000.0;

    // ---------------------- Further Task 1 -------------------------
    Nx--;
    //T1_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy);

    // ---------------------- Further Task 2 -------------------------
    // need to add choice of static or dynamic
    double T = 1;
    double Nt = 5;
    double rho = 7850;

    T2_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);

    std::cout << "\n\nEnd of program" << std::endl;
    return 0;
}


