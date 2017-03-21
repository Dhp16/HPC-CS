# include <iostream>
# include <string>
# include <iomanip> 
# include <cmath>



int main()
{
    double* X = new double[10];
    for(int i = 0; i < 10; i++){
        X[i] = 5;
    }   
    std::cout << "\nX" << std::endl;
    for(int i = 0; i < 10; i++)
        std::cout << X[i]<< std::endl;

    double* Y = new double[10];
    for(int i = 0; i < 10; i++)
        Y[i] = 2;

    std::cout << "\nY" << std::endl;
    for(int i = 0; i < 10; i++)
        std::cout << Y[i]<< std::endl;

    for(int i = 0; i < 10; i++)
        Y[i] = X[i]; 

    std::cout << "\nY (assigned to X)" << std::endl;
    for(int i = 0; i < 10; i++)
        std::cout << Y[i]<< std::endl;


    for(int i = 0; i < 10; i++)
        X[i] = 22;
    
   
    std::cout << "\nY (X re-assigned)" << std::endl;
    for(int i = 0; i < 10; i++)
        std::cout << Y[i]<< std::endl;


    return 0;
}