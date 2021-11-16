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
    
    double dmin=0, dmax=5.;
    //int n_point = 500;
    double h1 = 0.1, h2=h1/2, h3=h1/4;
    
    
    std::vector< std::vector<double> > fields_vect1;
    std::vector< std::vector<double> > fields_vect2;
    std::vector< std::vector<double> > fields_vect3;
    //std::vector< std::vector<double> > fields_vect1a;
    //std::vector< std::vector<double> > fields_vect2a;
    //std::vector< std::vector<double> > fields_vect3a;
    //std::vector< std::vector<double> > diff1;
    //std::vector< std::vector<double> > diff2;
    //std::vector< std::vector<double> > diff3;
    
    //ghost points at the boundary
    int gl = 1;
    int gr = 1;
    
    // initial condition //
    
    std::vector<double(*)(double, double)> initial_conditions;
    initial_conditions.push_back(&initial_gauss_PI);
    initial_conditions.push_back(&initial_null);
    //initial_conditions.push_back(&initial_null);
    
    
    // 2.35
    std::vector<double> parameters_ic_vector;
    for(double i=1.34;i<1.36;i=i+0.01)
    {
        parameters_ic_vector.push_back(i);
    }
    
    
    
    initialize_fields(fields_vect1,dmin,dmax,h1,initial_conditions,parameters_ic_vector[0]);
    initialize_fields(fields_vect2,dmin,dmax,h2,initial_conditions,parameters_ic_vector[0]);
    initialize_fields(fields_vect3,dmin,dmax,h3,initial_conditions,parameters_ic_vector[0]);
    //initialize_fields(fields_vect1a,dmin,dmax,h1,initial_conditions,parameters_ic_vector[1]);
    //initialize_fields(fields_vect2a,dmin,dmax,h2,initial_conditions,parameters_ic_vector[1]);
    //initialize_fields(fields_vect3a,dmin,dmax,h3,initial_conditions,parameters_ic_vector[1]);
   
    
    
    // setup of the diffential operator functions of the specific differential equation
    std::vector< double(*)(int ,int ,std::vector<std::vector<double>> &,double ,double ,std::vector<double> &, double ,std::vector<double (*)(std::vector<double>,int,double)> &,double (*)(double,int,std::vector<std::vector<double>> &,int,int,double,double),double ,int ,double, int ) > R_vector;
    //R_vector.push_back(&advection_eq_right_going);
    R_vector.push_back(&wave_eq_compactified_PI);
    R_vector.push_back(&wave_eq_compactified_PHI);
    //R_vector.push_back(&wave_eq_compactified_phi);
    
   
    // setup of the boundary conditions 
    std::vector< boundary_conditions_function> b_func;
    b_func.push_back(&no_boundary_conditions_PI_hyp);
    b_func.push_back(&no_boundary_conditions_PHI_hyp);
    //b_func.push_back(&no_boundary_conditions_phi_hyp);
    
    
    
    // setting the time integration parameters
    double  dt1 = h1*0.4, dt2=h2*0.4, dt3=h3*0.4, integration_interval = dt1*200, step_to_save = 200 ;
    //cout<<"dt1 = "<<dt1<<"\ndt2 = "<<dt2<<"\ndt3 = "<<dt3<<endl; //print the time steps
    // parameters vector that may be required from the system
    double v = 1.;
    double A = 1.;
    double s = dmax;
    std::vector<double> parameters;
    parameters.push_back(s);
    
    
    
    /* output writing version 1.0
    // we open an output file in order to write the name of the output files
    string file_path = "log_comp/data6/";
    ofstream names_file;
    names_file.open("./data/"+file_path+"names_file.csv");
    names_file <<"name\n";
    names_file.close();
    */
    
    // writing the output in a file
    string file_path = "./data/log_comp/data15/";
    
    
    
     //single run output writing
    
    string name_file = "hyp_slice_ampl_"+to_string(parameters_ic_vector[0])+"_dx_"+to_string(h1)+".csv";
    ofstream myfile2;
    myfile2.open (file_path+"name_of_file");
    myfile2<<"names\n"<<file_path+name_file<<"\n";
    myfile2.close();

    file_path.append(name_file);
    ofstream myfile;
    myfile.open (file_path);
    myfile.close();
    
    
    
    // dissipation coefficient for the art. diss.
    double epsilon1 = 0.0;
    
    
    std::vector<double (*)(std::vector<double>,int,double)> Dx;
    Dx.push_back(first_der_second_order_centered);
    //Dx.push_back(first_der_second_order_forward);
    //Dx.push_back(first_der_second_order_backward);

    //Dx.push_back(first_der_fourth_order_centered);
    /*
    // --------- EVOLUTION OF THE FUNCTION --------- //
    for(int l=0;l<parameters_ic_vector.size();l++)
    {
        string file_path = "./data/log_comp/data14/";
        string name_file = "ampl_"+to_string(parameters_ic_vector[l])+"_dx_"+to_string(h1)+".csv";
        file_path.append(name_file);
        ofstream myfile;
        myfile.open (file_path);
        myfile.close();
        cout<<"amplitutde "<<parameters_ic_vector[l]<<endl;
        std::vector< std::vector<double> > fields_vect1;
        std::vector< std::vector<double> > fields_vect2;
        std::vector< std::vector<double> > fields_vect3;
        
        initialize_fields(fields_vect1,dmin,dmax,h1,initial_conditions,parameters_ic_vector[l]);
        initialize_fields(fields_vect2,dmin,dmax,h2,initial_conditions,parameters_ic_vector[l]);
        initialize_fields(fields_vect3,dmin,dmax,h3,initial_conditions,parameters_ic_vector[l]);
        
        MOL_RK4(fields_vect1,&onestep_RK4_1,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
        MOL_RK4(fields_vect2,&onestep_RK4_1,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
        MOL_RK4(fields_vect3,&onestep_RK4_1,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    }
    
    */
    MOL_RK4(fields_vect1,&onestep_RK4_1,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    MOL_RK4(fields_vect2,&onestep_RK4_1,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    MOL_RK4(fields_vect3,&onestep_RK4_1,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    
    /*
    MOL_RK4(fields_vect1a,&onestep_RK4_1,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    //MOL_RK4(fields_vect2a,&onestep_RK4_1,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
    //MOL_RK4(fields_vect3a,&onestep_RK4_1,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,Dx,file_path);
      
    
    multiple_parameters_run(parameters_ic_vector,initial_conditions,initialize_fields,dmin,dmax,h1,h2,h3,dt1,dt2,dt3,integration_interval,step_to_save,Dx,R_vector,b_func,parameters,onestep_RK4_1,gl,gr,ghost_point_extrapolation_4_ord_spherical_symmetry,artificial_dissipation_2_Husa,epsilon1,print_f,file_path,MOL_RK4);
    */
    
    
    
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
