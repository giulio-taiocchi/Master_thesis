// Epa! //
#include <iostream>
#include <fstream>
#include<vector>
#include <math.h>
#include "myfunctions.h"
#include "spline.h"
using namespace std;

int main() 
{
    // integration in time //
    
    // setup of the fields and their initial condition in an domain=[dmin,dmax], after a h1-discretization
    
    double dmin=0, dmax=5;
    //int n_point = 500;
    double h1 = 0.01, h2=h1/2, h3=h1/4;
    std::vector< std::vector<double> > fields_vect1;
    std::vector< std::vector<double> > fields_vect2;
    std::vector< std::vector<double> > fields_vect3;
    std::vector< std::vector<double> > diff1;
    std::vector< std::vector<double> > diff2;
    std::vector< std::vector<double> > diff3;

    
    // initial condition //
    
    std::vector<double(*)(double, std::vector<double> &)> initial_conditions;
    initial_conditions.push_back(&initial_gauss_PI);
    initial_conditions.push_back(&initial_null);
    initial_conditions.push_back(&initial_null);
    
    
    std::vector<double> parameters_ic;
    parameters_ic.push_back(10);
    
    initialize_fields(fields_vect1,dmin,dmax,h1,initial_conditions,parameters_ic);
    initialize_fields(fields_vect2,dmin,dmax,h2,initial_conditions,parameters_ic);
    initialize_fields(fields_vect3,dmin,dmax,h3,initial_conditions,parameters_ic);
    
    
    
    // setup of the diffential operator functions of the specific differential equation
    std::vector< double(*)(int ,int ,std::vector<std::vector<double>> &,double ,double ,std::vector<double> &, double ,std::vector<double (*)(std::vector<double>,int,double)> &,double (*)(double,int,std::vector<std::vector<double>> &,int,int,double,double),double ,int ,double, int ) > R_vector;
    //R_vector.push_back(&advection_eq_right_going);
    R_vector.push_back(&wave_eq_spherical_PI);
    R_vector.push_back(&wave_eq_PHI);
    R_vector.push_back(&wave_eq_function);
    
   
    // setup of the boundary conditions 
    std::vector< void(*)(std::vector<std::vector<double>> &,std::vector<std::vector<double>> &,double ,double , double , int ,int , int ,double,double ,std::vector<double (*)(std::vector<double>,int,double)> &,double (*)(double,int,std::vector<std::vector<double>> &,int,int,double,double),double ,int)> b_func;
    //b_func.push_back(&adv_boundaries_left);
    b_func.push_back(&radiative_outer_boundaries_PI_we);
    b_func.push_back(&no_boundary_conditions_PHI);
    b_func.push_back(&no_boundary_conditions_phi);
    
    
    
    // setting the time integration parameters
    double  dt1 = h1*0.4, dt2=h2*0.4, dt3=h3*0.4, integration_interval = dt1*2500, step_to_save = 250 ;
    //cout<<"dt1 = "<<dt1<<"\ndt2 = "<<dt2<<"\ndt3 = "<<dt3<<endl; //print the time steps
    // parameters vector that may be required from the system
    double v = 1;
    std::vector<double> parameters;
    parameters.push_back(v);
    
    //ghost points at the boundary
    int gl = 1;
    int gr = 1;
    
    // we open an output file in order to write the name of the output files
    string file_path = "wave_equation/data34/";
    ofstream names_file;
    names_file.open("./data/"+file_path+"names_file.csv");
    names_file <<"name\n";
    names_file.close();
    
    
    
    // dissipation coefficient for the art. diss.
    double epsilon1 = 0.0;
    
    
    std::vector<double (*)(std::vector<double>,int,double)> Dx;
    Dx.push_back(first_der_second_order_centered);
    //Dx.push_back(first_der_fourth_order_centered);
    
    // --------- EVOLUTION OF THE FUNCTION --------- //
    MOL_RK4(fields_vect1,&onestep_RK4_1,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    MOL_RK4(fields_vect2,&onestep_RK4_1,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    MOL_RK4(fields_vect3,&onestep_RK4_1,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    
    
    // Convergence tests are performed 
            
    double self_test = self_conv_test(fields_vect1[0],fields_vect2[0],fields_vect3[0], h1, h2, &norm,&norm_of_diff );
    //double conv = conv_test(fields_vect1[0],fields_vect2[0],&line,&init_func,h1,h2, &norm, &norm_of_diff,dmin,dmax,h1,integration_interval);
    cout<<"epsilon:"<<epsilon1<<endl;
    
    /*
    // we open an output file in order to write the names of the "diff" files
    ofstream names_file_diff;
    names_file_diff.open("./data/names_file_diff.csv");
    names_file_diff <<"name_diff\n";
    names_file_diff <<"./data/h1_"+to_string(h1)+"_diff.csv\n";
    names_file_diff <<"./data/h1_"+to_string(h2)+"_diff.csv\n";
    names_file_diff <<"./data/h1_"+to_string(h3)+"_diff.csv\n";
    names_file_diff.close();
    
    std::vector<double> theo1;
    init_func(theo1,dmin,dmax,h1,line,integration_interval);
    std::vector<double> diff1_1;
    diff_vector(diff1_1, fields_vect1[0],theo1);
    diff1.push_back(diff1_1);
    
    print_f(diff1,dmin,h1,"diff.csv","./data/h1_"+to_string(h1)+"_" );
    names_file_diff <<"./data/h1_"+to_string(h1)+"_diff.csv\n";
    
    std::vector<double> theo2;
    init_func(theo2,dmin,dmax,h2,line,integration_interval);
    std::vector<double> diff2_2;
    diff_vector(diff2_2, fields_vect2[0],theo2);
    diff2.push_back(diff2_2);
    
    print_f(diff2,dmin,h2,"diff.csv","./data/h1_"+to_string(h2)+"_");   
    
    std::vector<double> theo3;
    init_func(theo3,dmin,dmax,h3,line,integration_interval);
    std::vector<double> diff3_3; 
    diff_vector(diff3_3, fields_vect3[0],theo3);
    diff3.push_back(diff3_3);
    
    print_f(diff3,dmin,h3,"diff.csv","./data/h1_"+to_string(h3)+"_");
    */
    return 0;
} 
