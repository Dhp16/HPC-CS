
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

// make sure there is a node under the concentrated force

void T1_inputs(){
    std::cout << "Initialising Elements" << std::endl;
    double L = Double_Value_Input("L");
    int Nx = Nx_Input();
    double A = Double_Value_Input("A");
    double I = Double_Value_Input("I");
    double E = Double_Value_Input("E");

    double l = L/Nx;    // number of elements
    return;
}



