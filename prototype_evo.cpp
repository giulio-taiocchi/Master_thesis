// Epa! //
#include <iostream>
#include <fstream>
#include<vector>
#include <math.h>
#include "myfunctions.h"
#include "spline.h"
#include <chrono>
using namespace std;
using namespace std::chrono;

int main() 
{
    auto start = high_resolution_clock::now();
    cout.precision(10);
    
    // declare the path where the input/output files are collected   
    string file_path = "./data/model3/data0/"; 
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
    
    std::vector<double(*)(double, double)> initial_conditions;
    initial_conditions.push_back(&initial_gauss);    //PI1
    initial_conditions.push_back(&initial_null);        //PHI1       
    initial_conditions.push_back(&initial_null);        //phi1
    initial_conditions.push_back(&initial_null);        //PI2
    initial_conditions.push_back(&initial_null);        //PHI2       
    initial_conditions.push_back(&initial_unity);       //phi2
    
    
    
    // setup of the diffential operator functions of the specific differential equation
    std::vector< evolution_function > R_vector;
    R_vector.push_back(&model3_PI1);
    R_vector.push_back(&model3_PHI1);
    R_vector.push_back(&model3_phi1);
    R_vector.push_back(&model3_PI2);
    R_vector.push_back(&model3_PHI2);
    R_vector.push_back(&model3_phi2);
    
   
    // setup of the boundary conditions 
    int number_of_bc = 6;
    std::vector< boundary_conditions_function> b_func(number_of_bc);
    b_func[0] = (&radiative_outer_boundaries_PI1_m3);
    b_func[1] = (&no_boundary_conditions_PHI1_m3);
    b_func[2] = (&no_boundary_conditions_phi1_m3);
    b_func[3] = (&radiative_outer_boundaries_PI2_m3);
    b_func[4] = (&no_boundary_conditions_PHI2_m3);
    b_func[5] = (&no_boundary_conditions_PHI2_m3);
    
    
    // setting the time integration parameters
    double  dt1 = h1*0.4, dt2=h2*0.4, dt3=h3*0.4;
    integration_interval=dt1*integration_interval;
    
    //cout<<"dt1 = "<<dt1<<"\ndt2 = "<<dt2<<"\ndt3 = "<<dt3<<endl; //print the time steps
   
    
    // setting the derivative operators vector
    std::vector<double (*)(std::vector<double>,int,double)> Dx;
    Dx.push_back(&first_der_second_order_centered);
    
    
    
    // multiple run setting
    
    // --------- EVOLUTION OF THE FUNCTION --------- //
    
    multiple_parameters_run(parameters_ic_vector,initial_conditions,initialize_fields,dmin,dmax,h1,h2,h3,dt1,dt2,dt3,integration_interval,step_to_save,Dx,R_vector,b_func,parameters,onestep_RK4_1,gl,gr,ghost_point_extrapolation_2_ord_TEM_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,print_f,file_path,MOL_RK4,ord);
    /*
    MOL_RK4(std::vector< std::vector<double> > fields_vect,one_step_function one_step, double dx, std::vector<double> param, double dt, double interval,    double dmin,    double dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path)
    
    multiple_layer(fields_vect, one_step, dx, param, dt, interval, dmin, dmax, R_vect, bc, step_to_save,  print_f, gl, gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path)
    */
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken: "<< duration.count()/1000000. << " seconds" << endl;
    
    return 0;
} 
