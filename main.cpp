# include <iostream>                
# include <iomanip>                        
# include <array>                       
# include <algorithm>                
# include "math.h"                  
# include <string>             
# include <boost/program_options.hpp>   
# include <vector>      
          
# include <mpi.h>  
         
# include "tools.h"     
# include "task1.h"      
# include "task2.h"        
# include "task3.h"    
# include "task4.h"     
# include "task5.h"  
 
//http://www.netlib.org/scalapack/explore-html/d0/d92/pdgbsv_8f_source.html

namespace po = boost::program_options;  
       
int main(int argc, char* argv[]){       
    std::cout << "\n\n-------------------------------------------------------------------------------- " << std::endl;
    std::cout << "---------------------------- Dominic Pickford 01272723 ------------------------- " << std::endl;  
    std::cout << "--------------------- High Performance Computing Assignment -------------------- " << std::endl;
    std::cout << "-------------------------------------------------------------------------------- " << std::endl;           

    const double b = 0.10;                        
    const double h = 0.120;           
     
    po::options_description desc("Calculates the dynamics of a Beam :D"); 
    desc.add_options()         
        ("L", po::value<double>(), "Length.")   
        ("Nx", po::value<int>(), "Number of elements.")    
        ("A", po::value<double>(), "Area")         
        ("I", po::value<double>(), "Moment of Inertia")     
        ("E", po::value<double>(), "Young's Modulus")          
                           
        ("T", po::value<double>(), "Total Time")           
        ("Nt", po::value<int>(), "Number of time steps")        
        ("rho", po::value<double>(), "Density") 
        ("equation", po::value<char>(), "Equation")               
        ("scheme", po::value<char>(), "Scheme")
        ("parallel", po::value<bool>(), "Run in Parallel"); 
        
    po::variables_map vm;          
    po::store(po::parse_command_line(argc, argv, desc), vm);     
    po::notify(vm); 
                  
    double  default_L = 10.0;     
    int     default_Nx = 24;   
    double  default_A = b*h;    
    double  default_I = b*pow(h, 3)/12; 
    double  default_E = 210000000000;   
            
    double  default_T = 55.0;  
    int     default_Nt = 200;         
    double  default_rho = 7850;    
    char default_equation = 'n';    // or dynamic
    char default_scheme = 'n';      // or implicit
    bool default_parallel = false;
       
    // read in all as strings and then validate inputs 
     
    // Inputs for Task 1     
    double L = vm.count("L")                  ?   vm["L"].as<double>()     : default_L;   
    int Nx = vm.count("Nx")                   ?   vm["Nx"].as<int>()       : default_Nx; 
    double A = vm.count("A")                  ?   vm["A"].as<double>()     : default_A; 
    double I = vm.count("I")                  ?   vm["I"].as<double>()     : default_I;  
    double E = vm.count("E")                  ?   vm["E"].as<double>()     : default_E;             
    // Inputs for Task 2           
    double T = vm.count("T")                  ?   vm["T"].as<double>()     : default_T;  
    int Nt = vm.count("Nt")                   ?   vm["Nt"].as<int>()       : default_Nt;  
    double rho = vm.count("rho")              ?   vm["rho"].as<double>()   : default_rho; 
    char equation = vm.count("equation")      ?   vm["equation"].as<char>() : default_equation;
    char scheme = vm.count("scheme")          ?   vm["scheme"].as<char>()   : default_scheme; 
    bool parallel = vm.count("parallel")      ?   vm["parallel"].as<bool>() : default_parallel;
    
    std::cout << "Scheme: " << scheme << std::endl;

    // ------------------ Task 1 (a): Inputs and Validation ----------
    if(L <= 0){
        std::cout << "NOTICE: A negative or null length is not appropriate, the default length of " << default_L << ", was selected." << std::endl;
        L = default_L;
    }
    if(Nx <= 0 || Nx % 2 != 0){
        std::cout << "NOTICE: A negative, null or odd number of elements is not appropriate, the default number of elements, " << default_Nx  << ", was selected." << std::endl;
        Nx = default_Nx;
    }
    if(A <= 0){
        std::cout << "NOTICE: A negative or null area is not appropriate, the default area, " << default_A  << ", was selected." << std::endl;
        A = default_A;
    }
    if(I <= 0){
        std::cout << "NOTICE: A negative or null moment of inertia is not appropriate, the default moment of inertia, " << default_I  << ", was selected." << std::endl;
        I = default_I;
    }
    if(E <= 0){
        std::cout << "NOTICE: A negative or null modulus of elasticity is not appropriate, the default moment of inertia, " << default_E  << ", was selected." << std::endl;
        E = default_E;
    }
    if(T <= 0){
        std::cout << "NOTICE: A negative or null total time is not appropriate, the default moment of inertia, " << default_T  << ", was selected." << std::endl;
        T = default_T;
    }
    if(Nt <= 0){
        std::cout << "NOTICE: A negative or null number of time steps is not appropriate, the default number of time steps, " << default_Nt << ", was selected." << std::endl;
        Nt = default_Nt;
    }
    if(rho <= 0){
        std::cout << "NOTICE: A negative or null density is not appropriate, the default density, " << default_rho  << ", was selected." << std::endl;
        rho = default_rho;
    }
    std::string Equation;
    if(equation == 's')
        Equation = "static";
    else if(equation == 'd')
        Equation = "dynamic";
    else
        Equation = "neither";
    std::string Scheme;
    if(scheme == 'e')
        Scheme = "explicit";
    else if(scheme == 'i')
        Scheme = "implicit";
    else 
        Scheme = "neither";
    if(Equation == "neither" && Scheme == "neither")
        std::cout << "NOTICE: The set parameters lead to none of the tasks being launched, please adjust your parameters for the scheme or the equation." << std::endl;

    // ----------------------- MPI Initialisation -------------------
    if(parallel){         
        int retval;     
        retval = MPI_Init(&argc, &argv);   
        if(retval == MPI_SUCCESS)   
           std::cout << "MPI Successfully Initialised" << std::endl;       // outputting twice, once for each process!
    }    

    // ----------------------- Variables needed ---------------------- 
    double l = L/Nx;    
    const double Qx = 0.0;       // no axial force in Task 1;      
    const double Qy = -1000.0;     //  originally in N/mm     
    const double Fy = -1000.0;                   
                        
    // ---------------------- Further Task 1 -------------------------
    Nx--;         
        
    std::cout <<"Equation: " << Equation << std::endl;
    std::cout <<"Scheme: " << Scheme << std::endl;
    std::cout <<"Parallel: " << parallel << std::endl;

    // ------ Launching each task depending on the parameters -------- 


    if(Equation == "static" && Scheme == "neither")   
        T1_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy);  
 
    if(Equation == "dynamic" && Scheme == "explicit" && parallel == false)   
        T2_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho); 
     
    if(Equation == "dynamic" && Scheme == "implicit" && parallel == false)  
        T3_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);      
 
    if(Equation == "dynamic" && Scheme == "explicit" && parallel == true)  
        T4_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);      
    
    if(Equation == "dynamic" && Scheme == "implicit" && parallel == true){
        Nx++;
        T5_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);
    }






    if(parallel) 
        MPI_Finalize();   
 
    std::cout << "\n\n-------------------------------------------------------------------------------- " << std::endl;   
    std::cout << "---------------------------- End of program ------------------------------------" << std::endl;
    std::cout << "--------------------------------------------------------------------------------" << std::endl;     
    return 0;    
} 
                            