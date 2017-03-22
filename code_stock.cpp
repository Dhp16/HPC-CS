    int Nx = Nx_Input();
    double A = Value_Input("A");
    double I = Value_Input("I");
    double E = Value_Input("E");

    int Nx_Input(){
    std::string input;
    int value(-1);
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

void central_Ke(const double A,const double E,const double I,const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[0] = 2*A*E/l;
    Ke[1+n] = 2*12*E*I/pow(l, 3);
    Ke[2+2*n] = 2*4*E*I/l;
    Ke[3+3*n] = Ke[0];
    Ke[4+4*n] = Ke[1+n];
    Ke[5+5*n] = Ke[2+2*n];

    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;

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
void first_Ke(const double A, const double E,const double I,const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[0] = A*E/l;
    Ke[1+n] = 12*E*I/pow(l, 3);
    Ke[2+2*n] = 4*E*I/l;
    //
    Ke[3+3*n] = 2*Ke[0];
    Ke[4+4*n] = 2*Ke[1+n];
    Ke[5+5*n] = 2*Ke[2+2*n];
    //
    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;
    //
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
void last_Ke(const double A, const double E,const double I, const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[3+3*n] = A*E/l;
    Ke[4+4*n] = 12*E*I/pow(l, 3);
    Ke[5+5*n] = 4*E*I/l;
    //
    Ke[0] = 2*Ke[3+3*n];
    Ke[1+1*n] = 2*Ke[4+4*n];
    Ke[2+2*n] = 2*Ke[5+5*n];
    //
    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;
    //
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


void central_Ke(const double A,const double E,const double I,const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[0] = 2*A*E/l;
    Ke[1+n] = 2*12*E*I/pow(l, 3);
    Ke[2+2*n] = 2*4*E*I/l;
    Ke[3+3*n] = Ke[0];
    Ke[4+4*n] = Ke[1+n];
    Ke[5+5*n] = Ke[2+2*n];

    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;

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
void first_Ke(const double A, const double E,const double I,const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[0] = A*E/l;
    Ke[1+n] = 12*E*I/pow(l, 3);
    Ke[2+2*n] = 4*E*I/l;
    //
    Ke[3+3*n] = 2*Ke[0];
    Ke[4+4*n] = 2*Ke[1+n];
    Ke[5+5*n] = 2*Ke[2+2*n];
    //
    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;
    //
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
void last_Ke(const double A, const double E,const double I, const double l){
    double Ke[36] = {0};
    int n = 6;
    Ke[3+3*n] = A*E/l;
    Ke[4+4*n] = 12*E*I/pow(l, 3);
    Ke[5+5*n] = 4*E*I/l;
    //
    Ke[0] = 2*Ke[3+3*n];
    Ke[1+1*n] = 2*Ke[4+4*n];
    Ke[2+2*n] = 2*Ke[5+5*n];
    //
    Ke[0+3*n] = -A*E/l;
    Ke[1+4*n] = -12*E*I/pow(l, 3);
    Ke[1+5*n] = 6*E*I/pow(l,2);
    Ke[2+4*n] = -Ke[1+5*n];
    Ke[2+5*n] = 2*E*I/l;
    //
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

void GFe_Take_Off_BCs(const double* G_Fe, double* G_FE_nBC, const int Nx){
    double* X = new double[Nx*6];

    for(int i = 0; i < Nx*6; i++)               // craeting a copy of G_Fe
        X[i] = G_Fe[i];

    back_array(X, 3, Nx);                       // shifting the values in the copy


    for(int i = 0; i < (Nx-1)*6; i++){          // assigning all the values in the copy to G_FE_nBC
        G_FE_nBC[i] = X[i];
    }

    delete X;
    return;
}

// -------------------------- Task 1 (c) ---------------------------------
void Global_Fe(double* G_Fe, double l, int Nx, double Qy){
    //double* G_Fe = new double[Nx*6];
    double Fe[6] = {0};
    pp_Fe(Fe, l, Nx, Qy, 1.);
    for(int i = 0; i < Nx;i++){
        G_Fe[i*6]   = Fe[0];
        G_Fe[i*6+1] = Fe[1];
        G_Fe[i*6+2] = Fe[2];
        G_Fe[i*6+3] = Fe[3];
        G_Fe[i*6+4] = Fe[4];
        G_Fe[i*6+5] = Fe[5];
    }
    return;
}
void banded_first_line(double* fl, const double* Ke, const double A, const double E,const double I, const double l, const int Nx){
    const int n = 6;
    //double* fl = new double[Nx*n];
    for(int i = 0; i < Nx*n; i++){
        if((i+1) % 6 == 0 && i != 0)
            fl[i] = Ke[1+5*n];      //6*E*I/pow(l,2);
        else
            fl[i] = 0; 
    }
    return;
}
void banded_second_line(double* sl, const double* Ke, const double A, const double E, const double I, const double l, const int Nx){
    const int n = 6;
    //double* sl = new double[Nx*n];
    sl[0] = 0;
    sl[1] = 0;
    sl[2] = 0;
    for(int i = 3; i < Nx*n; i = i+3){
        sl[i] = Ke[3*n];            // -A*E/l;
        sl[i+1] = Ke[1+4*n];        // -12*E*I/pow(l, 3);
        sl[i+2] = Ke[2+5*n];        // 2*E*I/l;
    }
    return;
}
void banded_third_line(double* th, const double* Ke, const double A, const double E, const double I, const double l, const int Nx){
    const int n = 6;
    th[5] = Ke[2+4*n];
    for(int i = 0; i < Nx*n; i++){
        if((i-4) % 6 == 0){
            th[i] = Ke[2+4*n];
        }
        else{
            th[i] = 0;
        }
    }     
}
void banded_fourth_line(double* fo, const double* Ke, const double A, const double E, const double I, const double l, const int Nx){
    const int n = 6;
    //double* fo = new double[Nx*n];
    for(int i = 0; i < Nx*n; i++){
        fo[i] = 0;
    }
    fo[2] = Ke[1+2*n];
    fo[Nx*n-1] = Ke[4+5*n];
}
void banded_mid_line(double* mid, const double* Ke, const double A, const double E, const double I, const double l, const int Nx){
    const int n = 6;
    //double* fi = new double[Nx*n];

    mid[0] = Ke[0];
    mid[1] = Ke[1+n];
    mid[2] = Ke[2+2*n];

    for(int i = 3; i < Nx*n-3; i = i+3){
        mid[i] = 2*Ke[0];
        mid[i+1] = 2*Ke[1+n];
        mid[i+2] = 2*Ke[2+2*n];
    }
    mid[Nx*n-3] = Ke[0];
    mid[Nx*n-2] = Ke[1+n];
    mid[Nx*n-1] = Ke[2+2*n];

    return;
}
void Global_Ke(double Ke[36], const int Nx, const double A, const double E, const double I, const double l){
    const int n = 6;

    double* G_Ke = new double[Nx*n*10];
    double* fl = new double[Nx*n];
    double* se = new double[Nx*n];
    double* th = new double[Nx*n];
    double* fo = new double[Nx*n];
    double* mid = new double[Nx*n];

    banded_first_line(fl, Ke, A, E, I, l, Nx);
    banded_second_line(se, Ke, A, E, I, l, Nx);
    banded_third_line(th, Ke, A, E, I, l, Nx);
    banded_fourth_line(fo, Ke, A, E, I, l, Nx);
    banded_mid_line(mid, Ke, A, E, I, l, Nx);

    int counter = 1;
    for(int i = 0; i < 9; i++){
        if(i == 0 || i == 8){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = fl[j];
            }
        }
        else if(i == 1 || i == 7){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = se[j];
            }
        }
        else if(i == 2 || i == 6){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = th[j];
            }
        }
        else if(i == 3 || i == 5){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = fo[j];
            }
        }
        else if(i == 4){
            for(int j = 0; j < Nx*n; j++){
                G_Ke[i*Nx*n+j] = mid[j];
            }
            back_array(fo, 1, Nx);
            back_array(th, 2, Nx);
            back_array(se, 3, Nx);
            back_array(fl, 4, Nx);
        }
    }

    for(int i = 0; i < 9; i++){
        for(int j = 0; j < Nx*n; j++){
            std::cout << std::setw(3);
            std::cout << G_Ke[i*Nx*n+j] << "  ";
        }
        std::cout << std::endl;
    }

    // for below diagonal put 0s at the end as opposed to the start

    return;
}
 


    for(int i = 0; i < Nx*n; i++){
        for(int j = 0; j < 13; j++){
            std::cout << std::setw(5);
            std::cout << N[i*Nx*n+j] << "  ";
        }
        std::cout << std::endl;
    }
    

        /*for(int i = 0; i < 13; i++){
        if(i == 8){
            for(int j = 0; j < Nx*n; j+=3){
                G_Mm[i*Nx*n + j] = 8;//2*1/2;
                G_Mm[i*Nx*n + j+1] = 8;//2*1/2;
                G_Mm[i*Nx*n + j+2] = 8;//2*alpha*pow(l, 2);
                counter = counter+3;
            }
        }
        else{
            for(int j = 0; j < Nx*n; j++)
                G_Mm[i*Nx*n + j] = 0; 
        }
    }*/
    for(int i = 0; i < 13*Nx*n; i++){
        if(i % 8 == 0){
            for(int j = 0; j < Nx*n; j++){
                G_Mm[i+j] = 8;//2*1/2;
                //G_Mm[i*Nx*n + j+1] = 8;//2*1/2;
                //G_Mm[i*Nx*n + j+2] = 8;//2*alpha*pow(l, 2);
                //counter = counter+3;
            }
        }
        else{
            //for(int j = 0; j < Nx*n; j++)
                G_Mm[i] = 0; 
        }
    }



void Global_Mass_Matrix_KS(const int Nx, const double rho, const double A, const double l){     // Global Mass Matrix K-Similar
    double alpha = 5; 
    int n = 3;
    double* G_Mm = new double[Nx*n*13];
    int counter = 0;

    /*std::cout << "Global Mass Matrix" << std::endl;
    for(int i = 0; i < Nx*n*13; i++){
       if(i % 13 == 0)
        std::cout << std::endl;
        std::cout << std::setw(5);
        std::cout << G_Mm[i] << " " ;
    }*/
    //std::cout << "COUNTER= " << counter << std::endl;

}

    for(int i = 0; i < Nx*n*9; i++){
        if(i == k*Nx*n){
            for(int i = 0; i < Nx*n){
                if(add == true)
                    AA[i] = AA[i] + M[counter]; 
                else
                    AA[i] = AA[i] - M[counter];
                counter++;
            }
        }
    }

        std::cout << "\nAfter" << std::endl;
    for(int i = 0; i < Nx*n*9; i++){
        if(i % (Nx*n) == 0)
            std::cout << "\n";
        std::cout << std::setw(15);
        std::cout << AA[i] << " ";
    }

        std::cout << "\nBefore" << std::endl;
    for(int i = 0; i < Nx*n*9; i++){
        if(i % 9 == 0)
            std::cout << "\n";
        std::cout << std::setw(15);
        std::cout << AA[i] << " ";
    }

        std::cout << "\nBefore" << std::endl;
    for(int i = 0; i < Nx*n*9; i++){
        if(i % (Nx*n) == 0)
            std::cout << std::endl;
        std::cout << output[n] << " ";
    }


    std::cout << "Global stiffness" << std::endl;
    for(int i = 0; i < 9+k; i++){
        for(int j = 0; j < Nx*n; j++){
            std::cout << std::setw(10);
            std::cout << G_Ke[i*Nx*n+j] << "  ";
        }
        std::cout << std::endl;
    }

        std::cout << "\nmid" << std::endl;
    for(int i = 0; i < Nx*3; i++)
        std::cout << bmid[i] << std::endl;



    // ------------------------ Optional --> Implement if finish  ------------------------------------
void bsfl(double* bsfl, const double* Ke, const int Nx){                         // banded first line
    int nK = 6;
    int n = 3;
    for(int i = 0; i < Nx*3; i++){
        if((i+1) % 2 == 0)
            bsfl[i] = Ke[1+nK];
        else
            bsfl[i] = 0;
    }
    std::cout << "\nbsfl" << std::endl;
    for(int i = 0; i < Nx*3; i++)
        std::cout << bsfl[i] << std::endl;
}
void bssl(double* bssl, const double* Ke, const int Nx){
    int nK = 6;
    int n = 3;
    for(int i = 0; i < Nx*3; i = i+3){
        bssl[i] = Ke[3*nK];
        bssl[i+1] = Ke[1+4*nK];
        bssl[i+2] = Ke[2+5*nK];
    }
    return;
}
void bstl(double* bstl, const double* Ke, const int Nx){
    int nK = 6;
    int n = 3;
    for(int i = 0; i < Nx*3; i++){
        if((i-1) % 3 == 0)
            bstl[i] = Ke[2+4*nK];
        else 
            bstl[i] = 0;
    }
}
void bsfol(double* bsfl,const double* Ke, const int Nx){
    for(int i = 0; i < Nx*3; i++)
        bsfl[i] = 0;
    return;
}
void bsmid(double* bsmid, const double* Ke, const int Nx){
    int nK = 6;
    for(int i = 0; i < Nx*3; i = i+3){
        bsmid[i] = 2*Ke[0];
        bsmid[i+1] = 2*Ke[1+nK];
        bsmid[i+2] = 2*Ke[2+2*nK]; 
    }
    return;
}

//////////////////////////////////////////////////////////////////////////////////
void BS_Global_Stiffness_Matrix(double* G_Ke, const double* Ke, const int Nx, const int k){         
    // B: Banded S: symetric    - Also Column Major
    const int n = 3;
    double* fl = new double[Nx*n];
    double* se = new double[Nx*n];
    double* th = new double[Nx*n];
    double* fo = new double[Nx*n];
    double* mid = new double[Nx*n];

    bsfl(fl, Ke, Nx);
    bssl(se, Ke, Nx);
    bstl(th, Ke, Nx);
    bsfol(fo, Ke, Nx);
    bsmid(mid, Ke, Nx);

    // make column major
    int counter = 0;
    for(int i = 0; i < Nx*n*5; i+=5){
        G_Ke[i] = fl[counter];
        G_Ke[i+1] = se[counter];
        G_Ke[i+2] = th[counter];
        G_Ke[i+3] = fo[counter];
        G_Ke[i+4] = mid[counter];
        counter++;
    }

    std::cout << "BandedSyemtic K GLOBAL" << std::endl;
    for(int i = 0; i < Nx*n*5; i++){
        if(i % 5 == 0)
            std::cout << "\n";
        std::cout << std::setw(10);
        std::cout << G_Ke[i] << "  ";
    }
}


    for(double i = 0; i < Nx/2; i++){
        double x = i/Nx*L;
        double u_distributed = Qy*pow(x, 2)*pow((L-x), 2)/(24*E*I);
        double u_concentrated = Fy*pow(x, 2)*(3*L-4*x)/(48*E*I);
        double u_total = u_distributed + u_concentrated;

        x_locs.push_back(x);
        asol.push_back(u_total);
    }
    int index = asol.size()-1;
    double x = L/2;
    double center = Qy*pow(x, 2)*pow((L-x), 2)/(24*E*I) + Fy*pow(x, 2)*(3*L-4*x)/(48*E*I);
    x_locs.push_back(0.5);
    asol.push_back(center);
    /*for(int i = Nx/2; i > 0; i--){
        double x = i/Nx*L;
        x_locs.push_back(x);
        double u_distributed = Qy*pow(x, 2)*pow((L-x), 2)/(24*E*I);
        double u_concentrated = Fy*pow(x, 2)*(3*L-4*x)/(48*E*I);
        double u_total = u_distributed + u_concentrated;
        asol.push_back(u_total);
    }
    int counter = 0;
    for(int i = Nx/2; i < Nx; i++){
        double x = (double)i/Nx*L;
        x_locs.push_back(x);
        asol.push_back(asol[index-counter]);
        counter++;
    }
    */

    int Nx = 4;
    double b = 0.10;
    double h = 0.120;
    double L = 10.0;
    double I = b*pow(h, 3)/12;
    double E = 210000000000;





    
    
    
    

void T2_inputsV3(const int Nx, const double b, const double h, const double L, const double A, const double I, const double E, const double l, const double Qx, const double Qy, const double Fy,
    const double T, const double Nt, const double rho){ 

    std::cout << "\nStarting Task 2" << std::endl;
    const int n = 3;
    double delta_t = T/Nt;
    int sizeA = 10;

    double Fe[6] = {0};
    pp_Fe(Fe, l, Nx, Qy, Qx);
    int k = 0;
    double Ke[36] = {0};
    pp_Ke(Ke, A, E, I, l);

    // COEF1= 1/delta_t^2 * [M]
    double* M_inv = new double[Nx*n];
    Global_Mass_Matrix(M_inv, Nx, rho, A, l);
    double u = 1/pow(delta_t, 2);
    Matrix_by_Scalar(M_inv, u, Nx);                     // multiply Global Mass Matrix by the scalar
    Inverse_Diagonal_Matrix(M_inv, Nx);                 // inverse the diagonal matrix

    double* U_i = new double[Nx*n];
    initialise_dynamic_array(U_i, Nx, n, 1);
    double* U_i_minus = new double[Nx*n];
    initialise_dynamic_array(U_i_minus, Nx, n, 1);
    double* U_i_plus = new double[Nx*n];
    initialise_dynamic_array(U_i_plus, Nx, n, 1);
    
    // ------------------  COEF2= [K] - 2/dealtaT^2 [M] ------------------------
    double* G_Me = new double[Nx*n];            
    Global_Mass_Matrix(G_Me, Nx, rho, A, l);            // generate M
    u = -2/pow(delta_t, 2);
    Matrix_by_Scalar(G_Me, u, Nx);                      // multiply M by u 
    double* BB = new double[Nx*n*9];
    Global_Stiffness_Matrix(BB, Ke, Nx);                // generate & populate K
    Matrix_Add_Diagonal(BB, G_Me, Nx, 4);               // add the two together
    delete G_Me;
    double* COEF2 = new double[Nx*n*9];
    output_array(BB, "Coef2", Nx*3, 9);
    Matrix_Transformer(COEF2, BB, Nx);                  //// Wrong 
    delete[] BB;
    // -------------------------------------------------------------------------
 
    // ------------------ COEF3= -1/2delta_t^2 [M] -----------------------------
    double* COEF3 = new double[Nx*n];
    Global_Mass_Matrix(COEF3, Nx, rho, A, l);           // correcting previous changes to G_Me
    u = -1/pow(delta_t, 2);
    Matrix_by_Scalar(COEF3, u, Nx);     
    // -------------------------------------------------------------------------

    int t_steps = (int)std::round(Nt);                  // need a fucntion to calculate the force every time step
    std::cout << "\nStarting loop" << std::endl;
    for(int i = 1; i < Nt+1; i++){
        
        // declaring variables to store the results of calculations
        double* Y2 = new double[Nx*n];      // COEF2*u{i}
        double* Y3 = new double[Nx*n];      // COEF2*u{i-1}
        double* RHS = new double[Nx*n];     

        // Calculate F
        double t_now = i*delta_t; 
        double* Fi = new double[Nx*n];
        Global_Force_Vector(Fi, Fe, Nx, k, Fy, t_now/T);    // populating the force vector and scaling the forces by t_now/T

        int lda = 9; // number of rows
        F77NAME(dgbmv) ('n', Nx*n, Nx*n, 4, 4, 1.0, COEF2, lda, U_i, 1, 0.0, Y2, 1);
        //F77NAME(dgbmv) ('n', Nx*n, Nx*n, 0, 0, 1.0, COEF3, 1, U_i_minus, 1, 0.0, Y3, 1);
        Diagonal_by_Vector(Y3, COEF3, U_i_minus, Nx);

        adding_three_arrays(RHS, Fi, Y2, Y3, Nx);

        // mutiplying by the inverse:
        //F77NAME(dgbmv) ('n', Nx*n, Nx*n, 0, 0, 1.0, M_inv, Nx*n, RHS, 1, 0.0, U_i_plus, 1);
        Diagonal_by_Vector(U_i_plus, M_inv, RHS, Nx);

        if(i < 5){
            std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
            std::cout << "\nY2: " << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << " " <<Y2[i] << std::endl;
            }
            std::cout << "\nY3: " << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << " " <<Y3[i] << std::endl;
            }
            std::cout << "\nRHS: " << std::endl; 
            for(int i = 0; i < Nx*n; i++){
                std::cout << " " <<RHS[i] << std::endl;
            }
            std::cout << "\nDisplacement" << std::endl;
            for(int i = 1; i < Nx*n; i+=3){
                std::cout << " " << U_i[i] << std::endl;
            }

        }

        // assigning the values for the following time step
        for(int i = 0; i < Nx*n; i++){
            U_i_minus[i] = U_i[i];
            U_i[i] = U_i_plus[i];
        }

        delete[] Fi;
        delete[] Y2;
        delete[] Y3;
        delete[] RHS;
    } 
    std::cout << "\nFinal Displacement at ts: "<< Nt << std::endl;
    for(int i = 1; i < Nx*n; i+=3)
            std::cout << " " << U_i[i] << std::endl;

    delete[] M_inv;
    delete[] COEF2;
    delete[] COEF3;
    delete[] U_i;
    delete[] U_i_plus;
    delete[] U_i_minus;
    return;
}



sub_from_vector(DIS_plus10, DIS10, vec_size);                                           // DIS_plus = DIS_plus - DIS_plus--> {u}n+1 = ({u}n+1 - {u}n)
        vector_by_scalar(DIS_plus10, 1/(beta*pow(delta_t, 2)), vec_size);                   // DIS_plus = DIS_plus*1/(beta*delta_t^2)       

        vector_by_scalar(VEL10, -1/(beta*delta_t), vec_size);
        vector_by_scalar(ACC10, 1-1/(2*beta), vec_size);


        // ----------------------- Equation 11 -------------------
        //double* ACC11 = new double[Nx*n];         Matrix_Copy(ACC11, ACC, vec_size);
        //double* ACC_plus11 = new double[Nx*n];  Matrix_Copy(ACC_plus11, ACC_plus, vec_size);

        for(int i = 0; i < vec_size; i++)
            VEL_plus[i] = VEL[i] + delta_t*(1-gamma)*ACC[i]+delta_t*gamma*ACC_plus;

        //vector_by_scalar(ACC11, delta_t*(1-gamma), vec_size);
        //vector_by_scalar(ACC_plus11, delta_t*gamma, vec_size);
        //for(int i = 0; i < vec_size; i++)
        //  ACC_plus[i] = VEL[i] + ACC11[i] + ACC_plus11[i];


                // ---------------------- Equation 10 -------------------
        double* DIS10 = new double[vec_size];       Matrix_Copy(DIS10, DIS, vec_size);
        double* DIS_plus10 = new double[vec_size];  Matrix_Copy(DIS_plus10, DIS_plus, vec_size);
        double* VEL10 = new double[vec_size];       Matrix_Copy(VEL10, VEL, vec_size);
        double* ACC10 = new double[vec_size];       Matrix_Copy(ACC10, ACC, vec_size);

        for(int i = 0; i < vec_size; i++)
            ACC_plus[i] = 1/(beta*pow(delta_t, 2))*(DIS_plus10[i] - DIS10[i]) - 1/(beta*delta_t)*VEL10[i] - (1/(2*beta))*ACC10[i];  /* ACC_plus = 1/(beta*delta_t^2)*({u}n+1 - {u}n) - 1/(beta*delta_t)*{udot}n - (1/2*beta -1){u dotdot}n*/
        
        delete[] DIS10;
        delete[] DIS_plus10;
        delete[] VEL10;
        delete[] ACC10; 

                // initialisation
        double* DIS12 = new double[vec_size]; Matrix_Copy(DIS12, DIS, vec_size);
        double* VEL12 = new double[vec_size]; Matrix_Copy(VEL12, VEL, vec_size);
        double* ACC12 = new double[vec_size]; Matrix_Copy(ACC12, ACC, vec_size);


        if(i == 1){                                                  
            std::cout << "\nK_effective" << std::endl;
            for(int i = 0; i < vec_size*13; i++){
                if(i % 13 == 0 && i != 0)
                    std::cout << "\n";
                std::cout << std::setw(12) << K_effective[i] << "   ";
            }
        }




    std::cout << "\nCOEF1" << std::endl;
    for(int i = 0; i < Nx*n; i++){
        std::cout << std::setw(10) << COEF1[i] << "         ";
    }

    std::cout << "\nCOEF3" << std::endl;
    for(int i = 0; i < Nx*n; i++){
        std::cout << std::setw(10) << COEF3[i] << "         ";
    }

    std::cout << "\nCOEF2" << std::endl;
    for(int i = 0; i < 9*Nx*n; i++){
        if(i % 9 ==0)
            std::cout << "\n";
        std::cout << std::setw(10) << COEF2[i] << "         ";
    }


    // Task 2 output
    if(i < 10 || i % 1000 == 0){
            std::cout << "\n\n---------  Time step " << i << "  ----------" << std::endl;
            std::cout << "T_now: " << t_now << std::endl;
            std::cout << " ----- Fi -----" << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << Fi[i] << " ";
            }
            std::cout << "\n ----- Y2 ----- " << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << Y2[i] << " ";
            }
            std::cout << "\n ----- Y3 ----- " << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << Y3[i] << " ";
            }
            std::cout << "\n ----- result ----- " << std::endl;
            for(int i = 0; i < Nx*n; i++){
                std::cout << result[i] << " ";
            }
            std::cout << "\n ----- U_plus ----- " << std::endl;
            for(int i = 1; i < Nx*n; i+=3){
                std::cout << U_plus[i] << " ";
            }        
        }


        if(i < 5){
            std::cout << "\n ----- Cmpy U_now going in ----- " << std::endl;
            for(int i = 1; i < Nx*n; i++){
                std::cout << copy_U_now[i] << " ";
            }   
        }

            std::cout << "\n\n\n------ Final Displacement TS: " << Nt << "------" << std::endl;
    for(int i = 1; i < Nx*n; i+=3){
        std::cout << U_plus[i] << " ";
    }       


    std::cout << "\nGlobal Mass Matrix" << std::endl; print(G_Mm, vec_size);

    std::cout << "\nK_effective" << std::endl;
    for(int i = 0; i < vec_size*13; i++){
        if(i % 13 == 0 && i != 0)
            std::cout << "\n";
        std::cout << std::setw(12) << K_effective[i] << "   ";
    }

        std::cout << "\nClean K" << std::endl;
    for(int i = 0; i < vec_size*13; i++){
        if(i % vec_size == 0 && i != 0)
            std::cout << "\n";
        std::cout << std::setw(12) << K_effective[i] << "   ";
    }

        // time step info in task 4
        if(rank == 1 && i < 15){
            std::cout << "\n-------------- Time step: " << i << " Tnow: " << t_now  << "-----------" << std::endl;
            std::cout << "Force" << std::endl;
            for(int i = 1; i < PARA_size; i+=3)
                std::cout << PARA_G_Fe[i] << "   ";
            std::cout << "\nU_Now: " << std::endl;
            for(int i = 1; i < PARA_size; i+=3)
                std::cout << U_now[i] << "   ";
        }

            std::cout <<"\nSOLUTION MESSAGE RECEIVED" << std::endl;
    for(int i = 0; i < PARA_size; i++)
        std::cout << Solution_Message[i] << "   ";

        std::cout << "\nPassed On: " << std::endl;
    std::cout << "\nVEC SIZE: " << std::endl;


     
  // if(rank == 1){
  //   std::cout << "\n\n\n\n\n\nRAnk: " << rank << " solution" << std::endl;
  //   for(int i = 1; i < PARA_size; i+=3)
  //     std::cout << U_now[i] << "   ";
  // }

  // combine solutions at the end   


     
    if(rank == 0)
        MPI_Send(Parcel, 6, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    else if(rank == 1){
        MPI_Send(Parcel, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
         
    if(rank == 0){ 
        MPI_Recv(&post_box, 6, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int counter = 0;
        for(int i = PARA_size-6; i < PARA_size; i++){
          U_now[i] = post_box[counter];
          counter++;
        }
    }
    else if(rank == 1){
        MPI_Recv(&post_box, 6, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i = 0; i < 6; i++)
          U_now[i] = post_box[i];
    } 
    if(rank == 0 && i == 9999){
      std::cout << "\nU_now After Change: " << std::endl;
      for(int i = 1; i < PARA_size; i+=3)
         std::cout << U_now[i] << "    ";
    }

    /*std::cout << "\nPARA_Ke final CM" << std::endl;
    for(int i = 0; i < PARA_nodes*n*(9+k); i++){
        if(i % (9+k) == 0)
            std::cout << std::endl;
        std::cout << std::setw(12) << PARA_Ke[i];  
    }            // Number of rows of 0s    */


        /*std::cout << "\nPARA_Ke with Banded 0s " << std::endl;
    for(int i = 0; i < (PARA_nodes+4)*n*(9+k); i++){
        if(i % ((PARA_nodes+4)*n) == 0)
            std::cout << std::endl;
        std::cout << std::setw(12) << PARA_Ke_W0s[i];  
    }*/