// Epa! //
#include <iostream>
#include <fstream>
#include<vector>
#include <math.h>
#include "myfunctions_parallel_ver.h"
#include "spline.h"
#include <chrono>
using namespace std;
using namespace std::chrono;
int main(int argc, char **argv) 
{
    auto start = high_resolution_clock::now();
    // set the decimal numbers precision of the output
    cout.precision(10);
    // calling the mpi initialization function
    int mynode, totalnodes;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    
    // declare the path where the input/output files are collected   
    string file_path = "./data/hyperboloidal_model3/data3/"; 
    // declare the name of the parameters file
    string parameter_file_name = "parameters_file_0";
    
    // initialize parameters
    double dmin,dmax,h1;
    
    // time interval of the integration and time steps to save in the output
    double integration_interval;
    int step_to_save;
    
    // number of left/right ghost points
    int gl,gr;
    
    // ord of the art. dissipation
    int ord;
    
    // dissipation coefficients for the art. diss.
    vector<double> epsilon1;
    
    // vector of parameters that enter in the initial conditions
    std::vector<double> parameters_ic_vector;
    
     // parameters vector that may be required from the system
    std::vector<double> parameters;
   
    // this function reads the parameters from a txt file
    read_parameters(file_path+parameter_file_name+".txt", dmin, dmax, h1, integration_interval, step_to_save, gl, gr, ord, epsilon1, parameters_ic_vector, parameters);

    
    // define the two higher resolutions tha are needed for the self convergence test
    double  h2=h1/2, h3=h1/4;
    
    // initial condition //
    
    std::vector<double(*)(double, double,vector<double>)> initial_conditions;
    initial_conditions_setting_hyp_Chi_m3_charvar_null_time_derivative(initial_conditions);
    
    
    
    // setup of the diffential operator functions of the specific differential equation
    std::vector< evolution_function > R_vector;
    evo_vector_setting_hyp_m3_Chi_charvar(R_vector);
    
    
    // setup of the boundary conditions 
    int number_of_bc = 6;
    std::vector< boundary_conditions_function> b_func(number_of_bc);
    bc_vector_setting_hyp_m3_Chi_charvar(b_func);
    
    
    // setting the time integration parameters
    double  dt1 = h1*0.4, dt2=h2*0.4, dt3=h3*0.4;
    integration_interval=dt1*integration_interval;
    
    //cout<<"dt1 = "<<dt1<<"\ndt2 = "<<dt2<<"\ndt3 = "<<dt3<<endl; //print the time steps
   
    
    // setting the derivative operators vector
    std::vector<double (*)(std::vector<double>,int,double)> Dx;
    Dx.push_back(&first_der_second_order_centered);
    
    
    
    
    // multiple run setting
    
    // --------- EVOLUTION OF THE FUNCTION --------- //
    
    multiple_parameters_run(parameters_ic_vector,initial_conditions,initialize_fields,dmin,dmax,h1,h2,h3,dt1,dt2,dt3,integration_interval,step_to_save,Dx,R_vector,b_func,parameters,onestep_RK4_1,gl,gr,ghost_point_extrapolation_2_ord_TEM_spherical_symmetry_charvar_Chi,artificial_dissipation_2_Husa,epsilon1,print_f,file_path,MOL_RK4,ord,status,totalnodes,mynode,request,communication );
    
   
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Time taken by processor "<<mynode<<" : "<< duration.count()/60. << " minutes" << endl;
    
    MPI_Finalize();
    
    return 0;
} 
