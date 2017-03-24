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
# include "task5b.h" 
      
//http://www.netlib.org/scalapack/explore-html/d0/d92/pdgbsv_8f_source.html
  
     
namespace po = boost::program_options;  
      
int main(int argc, char* argv[]){      
    bool task1 = false; 
    bool task2 = false;    
    bool task3 = false;    
    bool task4 = false;         
    bool task5 = true;       
                  
    //--------------------------------------------------
    if(task4 || task5){         
        int retval;    
        retval = MPI_Init(&argc, &argv);   
        if(retval == MPI_SUCCESS)   
           std::cout << "_Hash(16#{0})_cPentagram_@Pentagon{}hacking/*initialised" << std::endl;       // outputting twice, once for each process!
    }       
    //--------------------------------------------------
   
    double b = 0.10;                       
    double h = 0.120;          
            
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
        ;//("Scheme", po::value<std::string>(), "Scheme");
                 
    po::variables_map vm;          
    po::store(po::parse_command_line(argc, argv, desc), vm);     
    po::notify(vm);      
                  
    double  default_L = 10.0;     
    int     default_Nx = 6; 
    double  default_A = b*h;  
    double  default_I = b*pow(h, 3)/12; 
    double  default_E = 210000000000;   
          
    double  default_T = 55.0;  
    int     default_Nt = 200;        
    double  default_rho = 7850;   
    std::string default_scheme = "explicit";    
     
    // read in all as strings and then validate inputs 
      
    // Inputs for Task 1     
    const double L = vm.count("L")                  ?   vm["L"].as<double>()     : default_L;   
    int Nx = vm.count("Nx")                         ?   vm["Nx"].as<int>()       : default_Nx; 
    double A = vm.count("A")                        ?   vm["A"].as<double>()     : default_A; 
    double I = vm.count("I")                        ?   vm["I"].as<double>()     : default_I;  
    const double E = vm.count("E")                  ?   vm["E"].as<double>()     : default_E;          
    // Inputs for Task 2       
    const double T = vm.count("T")                  ?   vm["T"].as<double>()     : default_T;  
    const int Nt = vm.count("Nt")                   ?   vm["Nt"].as<int>()       : default_Nt; 
    const double rho = vm.count("rho")              ?   vm["rho"].as<double>()   : default_rho; 
    // Inputs for Task 3  
    //std::string scheme = vm.count("scheme")         ?   vm["scheme"].as<std::string>(): default_scheme;
              
              
    // --------------------------------------------------------------- 
    
    // ---------------------------------------------------------------  
   
	// ----------------------- Task 1 (a) ---------------------------- 
    // need to validate the inputs  
    double l = L/Nx;     
    // ----------------------- Variables needed ---------------------- 
    const double Qx = 0.0;       // no axial force in Task 1;      
    const double Qy = -1000.0;     //  originally in N/mm    
    const double Fy = -1000.0;                   
                
    // ---------------------- Further Task 1 -------------------------
    Nx--;        
       
    if(task1)   
        T1_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy);  
 
    if(task2)   
        T2_inputsV2(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho); 
     
    if(task3)  
        T3_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);      
 
    if(task4)  
        T4_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);      
    
    if(task5){
        Nx++;
        T5_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);
    }

    if(task4 || task5) 
        MPI_Finalize();   
 
      
    std::cout << "\n\nEnd of program" << std::endl;  
    return 0;    
} 
                          