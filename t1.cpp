void output_array(double A[6][6]){
    // improve to take any 2D array
    for(int i = 0; i < 6; i++){
        for(int j = 0; j < 6; j++){
            std::cout << std::setw(10);
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
}

void Ke_generator(double A, double E, double I, double l){  // Populating the element stiffness matrix
    double Ke[6][6] = {0};
    Ke[0][0] = A*E/l;
    Ke[1][1] = 12*E*I/pow(l, 3);
    Ke[2][2] = 4*E*I/l;
    Ke[3][3] = A*E/l;
    Ke[4][4] = Ke[1][1];
    Ke[5][5] = Ke[2][2];
    Ke[0][3] = -A*E/l;
    Ke[1][2] = 6*E*I/pow(l,2);
    Ke[1][4] = -12*E*I/pow(l, 3);
    Ke[1][5] = 6*E*I/pow(l,2);
    Ke[2][4] = -Ke[1][5];
    Ke[2][5] = 2*E*I/l;
    Ke[4][5] = -Ke[1][5];

    for(int i = 0; i < 6; i++){                         // Filling reamining fields knowing the matrix is symetric
        for(int j = 0; j < 6; j++){
            if(i == j)
                continue;
            if(i > j)
                Ke[i][j] = Ke[j][i];
        }
    }
    output_array(Ke);
    return;
}


// Qx and Qy magnitudes of the axial and transverse distributed loads
void Fe_generator(double l, double Qx, double Qy){      // Populating the Element Force Vector
    double[6] Fe;
    Fe[0] = Qx/2;
    Fe[1] = Qy/2;
    Fe[2] = Qy*l/2;
    Fe[3] = Qx/2;
    Fe[4] = Qy/2;
    Fe[5] = -Qy*l/12;
    for(int i = 0; i < Fe.size(); i++){
        Fe[i] = Fe[i]*l;
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
        value = std::stoi(input);
        if (value <= 0.0) {
            std::cout << "Nx must be positive!" << std::endl;
        }
        else if(value % 2 == 0.){
            std::cout << "Nx must be impair." << std::endl;
        }
     } while (value <= 0.0 || value % 2 == 0.);
    return value;
}

void T1_inputs(){
    std::cout << "Initialising Elements" << std::endl;
    double L = Double_Value_Input("L");
    int Nx = Nx_Input();
    double A = Double_Value_Input("A");
    double I = Double_Value_Input("I");
    double E = Double_Value_Input("E");
    double l = L/Nx;    // number of elements

    Ke_generator(A, E, I, l);
    return;
}
