#include "t1.cpp"

// Tools 
void output_array(double A[6][6]);
double Double_Value_Input(std::string name);
int Nx_Input();

void T1_inputs();

// make sure there is a node under the concentrated force

void Fe_generator(double l, double Qx, double Qy);
void Ke_generator(double A, double E, double I, double l);

void Global_Stiffness_Matrix_assembler(){	// needs to take Ke as argument



}

