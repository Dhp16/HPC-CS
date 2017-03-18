# include "t1.h"
/*#include <string>

    // improve to take any 2D array
void output_array(double *A, std::string s, int size){
    std::cout << "Printing: " << s << std::endl;
    for(int i = 0; i < size*size; i++){
        if(i != 0 && i % size == 0)
            std::cout << "\n";

        std::cout << std::setw(10);
        std::cout << A[i] << " ";
    }
    std::cout << "\n";
}

// -------------------------- T1b -------------------------- \\
void populate_Fe(double *Fe, int i, double l, int Nx, double Qy, double Qx){      // Will be different for each element, Qx --> function of x
    Fe[0] = Qx/2;
    Fe[1] = Qy/2;
    Fe[2] = Qy*l/2;
    Fe[3] = Qx/2;
    Fe[4] = Qy/2;
    Fe[5] = -Qy*l/12;
    for(int i = 0; i < 6; i++){                                     // could create an operator for this
        Fe[i] = Fe[i]*l;
    }
    return; 
}

// -------------------------- T1c -------------------------- \\
// Global Force Vector
void Global_Fe(double l, int Nx, double Qy){
    double* G_Fe = new double[Nx*6];
    for(int i = 0; i < Nx;i++){
        double Fe[6] = {0};
        populate_Fe(Fe, i, l, Nx, Qy, 1);                 // need to pass through information + Qx
        G_Fe[i*6]   = Fe[0];
        G_Fe[i*6+1] = Fe[1];
        G_Fe[i*6+2] = Fe[2];
        G_Fe[i*6+3] = Fe[3];
        G_Fe[i*6+4] = Fe[4];
        G_Fe[i*6+5] = Fe[5];
    }
    return;
}


double Double_Value_Input(std::string name){
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
    int value;
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

void Analytical_Solution(){

}
void T1_inputs(){


    // ----------------------- T1a -----------------------------------\\
    std::cout << "Initialising Elements" << std::endl;
    double L = Double_Value_Input("L");
    int Nx = Nx_Input();
    double A = Double_Value_Input("A");
    double I = Double_Value_Input("I");
    double E = Double_Value_Input("E");
    double l = L/Nx;    // number of elements

    // ----------------------- Variables needed ----------------------\\
    double Qx = 10;     
    double Qy = 10^3;    //  originally in N/mm   
    double x = 1;

    // ----------------------- T1b -----------------------------------\\
    std::cout << "HERE!" << std::endl;
    Global_Fe(l, Nx, Qy);

    return;
}
*/

void T1_inputs(){
    std::cout << "Starting Task1" << std::endl;
    double L = Value_Input("L");
    int Nx = Nx_Input();
    double A = Value_Input("A");
    double I = Value_Input("I");
    double E = Value_Input("E");
    double l = L/Nx;    // number of elements

    // ----------------------- Variables needed ----------------------\\
    double Qx = 10;     
    double Qy = 10^3;    //  originally in N/mm   
    double x = 1;

    // ----------------------- T1b -----------------------------------\\
    //std::cout << "HERE!" << std::endl;
    //Global_Fe(l, Nx, Qy);

    return;
};

std::cout << "Dimensions: " << std::endl;
    std::cout << "Force: " << 3*Nx << "x" << 1 << std::endl;
    std::cout << "Stiffness: " << 3*Nx << "x" << 10 << std::endl;