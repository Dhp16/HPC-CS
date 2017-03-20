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
  
namespace po = boost::program_options;  
      
int main(int argc, char* argv[]){   
    bool task1 = false;
    bool task2 = false;  
    bool task3 = false;   
    bool task4 = true;      
            
    //--------------------------------------------------
    if(task4){       
        int retval;    
        retval = MPI_Init(&argc, &argv); 
        if(retval == MPI_SUCCESS)   
           std::cout << "MPI_SUCCESS!" << std::endl;       // outputting twice, once for each process!
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
    /*int rank; 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int message;

    if(rank == 0)
        message = 666;
    else  
        message = 333;
    int post_box; 

    if(rank == 0)
        MPI_Send(&message, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    else if(rank == 1){
        MPI_Recv(&post_box, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "Process 1 received the following number: " <<  post_box << std::endl;
    }
    */   
 
    // Successfully passing messages 
    /*
    double Info[12];
    if(rank == 0)
        for(int i = 0; i < 12; i++) Info[i] = 7777;
    
    else 
        for(int i = 0; i < 12; i++) Info[i] = 1111;
    
    double post_box[12]; 

    if(rank == 0)
        MPI_Send(&Info, 12, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    else if(rank == 1){
        MPI_Send(&Info, 12, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }


    if(rank == 0){
        MPI_Recv(&post_box, 12, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "Process 0 received the following numbers: " << std::endl;
        for(int i = 0; i < 12; i++)
            std::cout << rank <<"   " <<  post_box[i] << "  ";
    }
    else if(rank == 1){
        MPI_Recv(&post_box, 12, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::cout << "Process 1 received the following numbers: " << std::endl;
        for(int i = 0; i < 12; i++)
            std::cout << post_box[i] << "  ";
    }
    */
    // ---------------------------------------------------------------
 
	// ----------------------- Task 1 (a) ----------------------------
    // need to validate the inputs 
    double l = L/Nx;    
    // ----------------------- Variables needed ----------------------
    const double Qx = 0.0;       // no axial force in Task 1;     
    const double Qy = -1000.0;    //  originally in N/mm   
    const double Fy = -1000.0;              
               
    // ---------------------- Further Task 1 -------------------------
    Nx--;      
  
    if(task1) 
        T1_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy);  
   
    // ---------------------- Further Task 2 -------------------------
    // need to add choice of static or dynamic  

    std::string scheme = "explicit";

    if(task2)  
        T2_inputsV2(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho); 
     
    if(task3)
        T3_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);
 
    if(task4) 
        T4_inputs(Nx, b, h, L, A, I, E, l, Qx, Qy, Fy, T, Nt, rho);
    
    if(task4) 
        MPI_Finalize();  
    
    std::cout << "\n\nEnd of program" << std::endl; 
    return 0;
} 
   
    
                                