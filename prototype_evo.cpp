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
    
    
    // setup of the fields and their initial condition in an domain=[dmin,dmax]
    
    double dmin=0, dmax=2.5;
    //int n_point = 500;
    double h1 = 0.05, h2=h1/2, h3=h1/4;
    
    //std::vector< std::vector<double> > diff1;
    //std::vector< std::vector<double> > diff2;
    //std::vector< std::vector<double> > diff3;
    
    //ghost points at the boundary
    int gl = 1;
    int gr = 1;
    
    // dissipation coefficient for the art. diss.
    //vector<double> epsilon1 = {0.045,0.044,0.02,0.008}; 
    vector<double> epsilon1;
    for (double e=0;e<0.02;e=e+0.01)
    {
        epsilon1.push_back(e);
    }
    int ord = 2;
    
    // initial condition //
    
    std::vector<double(*)(double, double)> initial_conditions;
    initial_conditions.push_back(&initial_gauss_PI);
    initial_conditions.push_back(&initial_null);
    //initial_conditions.push_back(&initial_null);
    
    
    // 2.35
    std::vector<double> parameters_ic_vector;
    for(double a=1;a<1.1;a=a+1)
    {
        parameters_ic_vector.push_back(a);
    }
    
    
    
    
    // setup of the diffential operator functions of the specific differential equation
    std::vector< evolution_function > R_vector;
    //R_vector.push_back(&advection_eq_right_going);
    R_vector.push_back(&model1_PI);
    R_vector.push_back(&wave_eq_PHI);
    //R_vector.push_back(&wave_eq_compactified_phi);
    
   
    // setup of the boundary conditions 
    int number_of_bc = 2;
    std::vector< boundary_conditions_function> b_func(number_of_bc);
    b_func[0] = (&radiative_outer_boundaries_PI_m1);
    b_func[1] = (&no_boundary_conditions_PHI);
    //b_func[1] = (&no_boundary_conditions_phi);
    
    
    
    // setting the time integration parameters
    double  dt1 = h1*0.4, dt2=h2*0.4, dt3=h3*0.4, integration_interval = dt1*1200, step_to_save = 300 ;
    //cout<<"dt1 = "<<dt1<<"\ndt2 = "<<dt2<<"\ndt3 = "<<dt3<<endl; //print the time steps
    // parameters vector that may be required from the system
    double v = 1.;
    double A = 1.;
    double s = dmax;
    std::vector<double> parameters;
    parameters.push_back(A);
    
    // setting the derivative operators vector
    std::vector<double (*)(std::vector<double>,int,double)> Dx;
    Dx.push_back(&first_der_second_order_centered);
    
    
    
    
    // writing the output in a file
    string file_path = "./data/log_comp/data16/";   
    
    
    //single run 
    /*
    // output files setting
    string name_file = to_string(parameters_ic_vector[0])+"_dx_"+to_string(h1)+".csv";
    ofstream myfile2;
    myfile2.open (file_path+"name_of_file");
    myfile2<<"names\n"<<file_path+name_file<<"\n";
    myfile2.close();

    file_path.append(name_file);
    ofstream myfile;
    myfile.open (file_path);
    myfile.close();
    
    // initialization of the fields
    std::vector< std::vector<double> > fields_vect1 = initialize_fields(dmin,dmax,h1,initial_conditions,parameters_ic_vector[0],gl,gr,ord);
    std::vector< std::vector<double> > fields_vect2 = initialize_fields(dmin,dmax,h2,initial_conditions,parameters_ic_vector[0],gl,gr,ord);
    std::vector< std::vector<double> > fields_vect3 = initialize_fields(dmin,dmax,h3,initial_conditions,parameters_ic_vector[0],gl,gr,ord);
    
    
    // evolution
    MOL_RK4(fields_vect1,&onestep_RK4_1,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,ord,Dx,file_path);
    MOL_RK4(fields_vect2,&onestep_RK4_1,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,ord,Dx,file_path);
    MOL_RK4(fields_vect3,&onestep_RK4_1,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,ord,Dx,file_path);
    */
    
    
    // multiple run setting
    
    // --------- EVOLUTION OF THE FUNCTION --------- //
    
    multiple_parameters_run(parameters_ic_vector,initial_conditions,initialize_fields,dmin,dmax,h1,h2,h3,dt1,dt2,dt3,integration_interval,step_to_save,Dx,R_vector,b_func,parameters,onestep_RK4_1,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,print_f,file_path,MOL_RK4,ord);
    
    
    // Convergence tests are performed 
            
    //double self_test = self_conv_test(fields_vect1[0],fields_vect2[0],fields_vect3[0], h1, h2, &norm,&norm_of_diff );
    //double conv = conv_test(fields_vect1[0],fields_vect2[0],&line,&init_func,h1,h2, &norm, &norm_of_diff,dmin,dmax,h1,integration_interval);
    //cout<<"epsilon:"<<epsilon1<<endl;
    
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
