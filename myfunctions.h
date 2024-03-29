#include <math.h>
#include <fstream>
#include <sstream>
#include "spline.h"
using namespace std;

typedef void (*print_function)(std::vector< std::vector<double> > fields_vect, double dmin, double dx, string name_file,string name_folder, int gl, int gr);


typedef double (*artificial_dissipation_function)(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt) ;

typedef std::vector<double (*)(std::vector<double>,int,double)>  derivative_vector;

typedef void(*ghost_point_extrapolation_function)(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax);

typedef double (*evolution_function)(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl);

typedef void(*boundary_conditions_function)(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param);

typedef vector<vector<double>>(*one_step_function)(std::vector< std::vector<double> > fields_vect,double dmin,double dmax,double dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord)   ;

typedef void(*method_of_line_function)(std::vector< std::vector<double> > fields_vect,one_step_function one_step, double dx, std::vector<double> param, double dt, double interval,    double dmin,    double dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path);

typedef vector<vector<double>>(*initialization_function)(double d_min,double d_max,double dx,std::vector<double(*)(double,double)> funcs,double param_ic,int gl, int gr,int ord);


void multiple_parameters_run(std::vector<double>& parameters_ic_vector, std::vector<double(*)(double, double)>& initial_conditions, initialization_function initialize_fields, double dmin, double dmax, double h1, double h2, double h3, double dt1, double dt2, double dt3, double integration_interval,int step_to_save, derivative_vector & Dx ,std::vector< evolution_function > & R_vector, std::vector< boundary_conditions_function > & b_func,  std::vector<double> & parameters, one_step_function one_step, int gl,int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2, vector<double> epsilon1,print_function print_f, string file_path, method_of_line_function MOL_RK4,int ord)
{
    ofstream myfile2;
    //myfile2.open (file_path+"name_of_file");
    //myfile2<<"names\n";
    //myfile2.close();
    
    for(int e=0;e<epsilon1.size();e++)
    {
        for (int l=0;l<parameters_ic_vector.size();l++)
        {
            cout<<"\nepsilon = "<<epsilon1[e]<<endl;
            cout<<"\namplitude = "<<parameters_ic_vector[l]<<endl;
            // initialize fields
            
            
            
            std::vector< std::vector<double> > fields_vect1 = initialize_fields(dmin,dmax,h1,initial_conditions,parameters_ic_vector[l],gl,gr,ord);     
            std::vector< std::vector<double> > fields_vect2 = initialize_fields(dmin,dmax,h2,initial_conditions,parameters_ic_vector[l],gl,gr,ord);
            std::vector< std::vector<double> > fields_vect3 = initialize_fields(dmin,dmax,h3,initial_conditions,parameters_ic_vector[l],gl,gr,ord);
            // setting the output variables
            //cout<<"\n lun vec \n"<<endl;
            string name_file = "ampl_"+to_string(parameters_ic_vector[l])+"_eps"+to_string(epsilon1[e])+"_dx_"+to_string(h1)+"steps"+to_string(step_to_save)+"last_time"+to_string(integration_interval)+".csv";
            
            string complete_path = file_path+name_file;
            ofstream myfile;
            
            myfile.open (complete_path);
            myfile.close();
            
            MOL_RK4(fields_vect1,one_step,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path);
            MOL_RK4(fields_vect2,one_step,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path);
            MOL_RK4(fields_vect3,one_step,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path);
        }
    }
}
    

void MOL_RK4(std::vector< std::vector<double> > fields_vect,one_step_function one_step, double dx, std::vector<double> param, double dt, double interval,    double dmin,    double dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path)
{   
    
    
    cout<<"--- Method of lines called ---\ndx = "<<dx<<"\ndt = "<<dt<<"\nDomain = ["<<dmin<<","<<dmax<<"]\nlast time objective :"<<interval<<"\n";
    
    if (ord>gl)
    {
        gl = ord;
    }
    if(ord>gr)
    {
        gr = ord;
    }
    
    
    int ns = interval / dt; // the RK will be called ns times
    int m = (ns/step_to_save); // every m step the data are printed to a file
    int counter = 0;
    //ofstream myfile_times;
    double last;
    for (double t=0;t<interval+dt/1.5;t=t+dt) 
    {        
        
        // we print to file only for some time
        if ((counter%m)==0)
        {
            
            //ofstream myfile;
            //myfile.open (file_path, ios::app);
            //myfile<<t<<"\n";
            //myfile.close();
            cout<<"print the fields at time"<<t<<endl;
            print_f(fields_vect,dmin,dx,file_path,file_path,gl,gr); // the print_f function is called
            //cout<<t<<endl;
        }
        fields_vect = one_step(fields_vect,dmin,dmax,dx,param,dt,R_vect,bc,t,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon,Dx,ord);
        
        counter += 1;
    last = t;
    }
    
    //cout<<"last time included:"<<last<<endl;
    
}





/*
void onestep_RK1_1(std::vector< std::vector<double> > &fields_vect,double dmin,double dmax,double dx,std::vector<double> &param, double dt, std::vector< evolution_function > &evo,std::vector< boundary_conditions_function > &bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx)  
{
    //cout<<"size: "<<fields_vect[0].size()<<endl;
    std::vector<double> support;
    int N = fields_vect.size();
     we copy the field vector, we will use it at the end of the function
    /*
    std::vector<std::vector<double>> copy_fields_vect;
    for (int j=0; j <N; j++)
    {
        copy_fields_vect.push_back(support);
        for (int i=0;i<fields_vect[j].size();i++)
        {
            copy_fields_vect[j].push_back(fields_vect[j][i]);
        }
    }  
    // A second order Kreis Oliger artificial dissipation is used
    
    
    if (ord>gl)
    {
        gl = ord;
    }
    if(ord>gr)
    {
        gr = ord;
    }
    
    
    
    
    std::vector<std::vector<double>> k1;
    /*
    // populate the ghost zone of the fields, we need them to calculate k1
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }
    
    // k1 building
    for (int j=0; j <N; j++)
    {
        k1.push_back(support);
        // evualuating the "bulk" of k1
        // we have to consider the physical domain of the fields (so exclude the GP) and then exclude the boundaries value (+1 and -1)
        for (int i=gl+1;i<(fields_vect[j].size()-1)-gr;i++)
        {
            
            k1[j].push_back(dt*evo[j](j,i,fields_vect,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl));
        }
        //cout<<"prima k1: "<<k1[j].size()<<endl;
        // evaluating the boundary of k1
        bc[j](k1,fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param);
        //cout<<"dopo k1: "<<k1[j].size()<<endl;


        //cout<<"size k1: "<<k1[j].size()<<"\nsize field vector: "<<fields_vect[j].size()<<endl;
        
        
    }
    
    
    // we create a new vector that contains all the new fields. It is a support vector that will be swapped with the old one
    
    std::vector< std::vector<double> > new_fields_vect;
        
    for (int j=0; j <N; j++)
    {
        new_fields_vect.push_back(support);
        for (int i=0;i<k1[j].size();i++)
        {
            new_fields_vect[j].push_back(fields_vect[j][i+gl] + k1[j][i+gl]);
        }
    }  
    
    
    //cout<<"old "<<fields_vect[0].size()<<"new "<<new_fields_vect[0].size()<<endl;
    for (int j=0; j <N; j++)
    {
        for (int i=0;i<k1[j].size();i++)
        {
            fields_vect[j][i]=new_fields_vect[j][i];
        }
    }
    //new_fields_vect.swap(fields_vect);
    
}*/

vector<vector<double>> onestep_RK4_1(std::vector< std::vector<double> > fields_vect,double dmin,double dmax,double dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord)  
{
    //cout<<"size: "<<fields_vect[0].size()<<endl;
    int N = fields_vect.size();
    //cout<<"number of fields:"<<N<<endl;
    int S =  int( ((dmax+dx*double(gr))-(dmin-dx*double(gl)) +dx/2 ) / dx) + 1;
    //cout<<"S: "<<S<<endl;
    // A second order Kreis Oliger artificial dissipation is used
    
    vector<vector<double>> k1(N, vector<double>(S)) ;
    std::vector<std::vector<double>> support_k1(N, vector<double>(S));
    std::vector<std::vector<double>> k2(N, vector<double>(S));
    std::vector<std::vector<double>> support_k2(N, vector<double>(S));
    std::vector<std::vector<double>> k3(N, vector<double>(S));
    std::vector<std::vector<double>> support_k3(N, vector<double>(S));
    std::vector<std::vector<double>> k4(N, vector<double>(S));
    std::vector<std::vector<double>> support_k4(N, vector<double>(S));
       
    // populate the ghost zone of the fields, we need them to calculate k1
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }
    // k1 building
    for (int j=0; j <N; j++)
    {
        // evualuating the "bulk" of k1
        // we have to consider the physical domain of the fields (so exclude the GP) and then exclude the boundaries value (+1 and -1)
        for (int i=gl+1;i<(fields_vect[j].size()-1)-gr;i++)
        {
            k1[j][i] = (evo[j](j,i,fields_vect,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl)  
                            +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
        }
        //cout<<"prima k1: "<<k1[j].size()<<endl;
        // evaluating the boundary of k1
        bc[j](k1,fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param);
        //cout<<"dopo k1: "<<k1[j].size()<<endl;

        ghost_point_extrapolation(k1, t,dx,dt,j,gl,gr,dmin,dmax);
        //cout<<"post extr k1:"<<k1[j].size()<<endl;

        //cout<<"k1"<< k1.size()<<endl;
    // computing the argument for the next coefficient ki
    
        for (int i=0;i<k1[j].size();i++)
        {
            support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];
        }
    }
    
    
    
    // k2 building
    for (int j=0; j <N; j++)
    {
        for (int i=1+gl;i<(support_k1[j].size()-1-gr);i++)
        {
            k2[j][i] = (evo[j](j,i,support_k1,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt,gl)+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
        }
        //cout<<"prima k2:"<<k2[j].size()<<endl;
        bc[j](k2,support_k1, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param);
        //cout<<"dopo k2:"<<k2[j].size()<<endl;
        
        ghost_point_extrapolation(k2, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax);
        //cout<<"post extr k2:"<<k2[j].size()<<endl;
        for (int i=0;i<k2[j].size();i++)
        {
            support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];
        }
    
    }
    /*
    // populate the ghost zone of support k2, we need them to calculate k3
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(support_k2, t+dt/2,dx,dt,j,gl,gr,dmin,dmax);
    }    
    */
    // k3 building
    for (int j=0; j <N; j++)
    {
        for (int i=gl+1;i<(k2[j].size()-1)-gr;i++)
        {
            k3[j][i] = (evo[j](j,i,support_k2,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt,gl)+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
        }
        bc[j](k3,support_k2, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param);
        
        
        ghost_point_extrapolation(k3, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax);
        
        for (int i=0;i<k3[j].size();i++)
        {
            support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];
        }
    }
   
   
    // k4 building    
    for (int j=0; j <N; j++)
    {
        for (int i=gl+1;i<(k3[j].size()-1-gr);i++)
        {
            k4[j][i] = evo[j](j,i,support_k3,dx,dmin,param,t+dt,Dx,artificial_diss,epsilon,ord,dt,gl)+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
        }
        bc[j](k4,support_k3, t+dt,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param);
        
        ghost_point_extrapolation(k4, t+dt,dx,dt,j,gl,gr,dmin,dmax);
    }
    // we create a new vector that contains all the new fields. It is a support vector that will be swapped with the old one
    
    std::vector<std::vector<double>> new_fields_vect(N, vector<double>(S));

    for (int j=0; j <N; j++)
    {
        for (int i=0;i<k4[j].size();i++)
        {
            new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.;
        }
    }    
    //cout<<"old "<<fields_vect.size()<<"new "<<new_fields_vect.size()<<endl;
    return(new_fields_vect);
}

/*
// runge kutta integrator for a system of equations as dy/dx=F
// where dy/dx and F are vectors
// input: F( vec(y_0,.....,y_n),x), y -> vector with initial values of y in different point (x), extreme of the integral [dmin,dmax], increment dx
void RK4(std::vector<double> y,std::vector<double(*)(std::vector<double>,double)> &F, double dmin, double dmax, double dx)
{
    
    std::vector<double> k1;
    std::vector<double> support_k1;
    std::vector<double> k2;
    std::vector<double> support_k2;
    std::vector<double> k3;
    std::vector<double> support_k3;
    std::vector<double> k4;
    std::vector<double> support_k4;
    
    int N = y.size();
    for (double x=dmin; x<=dmax; x=x+dx )
    {
        for (int i=0;i<N;i++)
        {
            k1.push_back(dx*F[i](y,x));
            support_k1.push_back(k1[i]/2+y[i]);
        }
        
        for (int i=0;i<N;i++)
        {
            k2.push_back(dx*F[i](support_k1,x+dx/2));
            support_k2.push_back(k2[i]/2+y[i],x+dx/2);
        }
        
        for (int i=0;i<N;i++)
        {
            k3.push_back(dx*F[i](support_k2,x+dx/2));
            support_k3.push_back(k3[i]+y[i]);
        }
        
        for (int i=0;i<N;i++)
        {
            k4.push_back(dx*F[i](support_k3,x+dx));
        }
        
        std::vector new_y;
        for (int i=0;i<N;i++)
        {
            new_y.push_back(y[i]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6);
        }
    }
    
}

*/

// ------------ TYPES OF D.E. --------------- //



// ----------- // ADVECTION EQUATION // ----------- //



// advection equation, central finite difference, 2nd order method
double advection_eq_left_going(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    {
        return (param[0]*Dx[0](fields_vect[0],ind_space,dx));
    }
}

double advection_eq_right_going(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    {
        return (param[0]*(-1.)*Dx[0](fields_vect[0],ind_space,dx));
    }
    
}

// ----------- // WAVE EQUATION // ----------- //


double wave_eq_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    return (Dx[0](fields_vect[1],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[1],ind_space,dx));

}
    
double wave_eq_PHI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[0],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double wave_eq_phi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[0][ind_space]);
}


//  wave equation in spherical symmetry
double wave_eq_spherical_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt,int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    // Dx[0](fields_vect[1],ind_space,dx) + 2./x * fields_vect[1][ind_space]
    // 3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) 
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) );
}

// R RESCALING HYPERBOLOIDAL foliation //
// here are reported the evolution equations for the auxiliary function of the rescaled (by R) and compactified wave equation
double wave_eq_compactified_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double rate_of_square =pow(r,2)/pow(s,2);
    //double R = r/(1-rate_of_square);
    //double Rprime = (1+rate_of_square)/pow((1-rate_of_square),2);
    //double Hprime = 1-1/Rprime;
    
    //w we introduce a variable for 1/(R'(1-H'^2))
    double coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variable for H'/(R'*(1-H'^2))
    double coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
    
    return (-coefficient2*Dx[0](fields_vect[0],ind_space,dx)
            -coefficient1*Dx[0](fields_vect[1],ind_space,dx)
            );
        
}

double wave_eq_compactified_PHI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double rate_of_square = pow(r,2)/pow(s,2);
    
    //w we introduce a variable for 1/(R'(1-H'^2))
    double coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variavle for H'/(R'*(1-H'^2))
    double coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
 
    return (-coefficient2*Dx[0](fields_vect[1],ind_space,dx)
            -coefficient1*Dx[0](fields_vect[0],ind_space,dx));
        
}

double wave_eq_compactified_phi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    return(-fields_vect[0][ind_space]);
}


//------------------- CHI RESCALING HYPERBOLOIDAL FOLIATION -------------------//

// compactification: T = t + H(R), H'=1-1/R', R=r/(1-(r/s)^2)
// rescaling: by chi(R) = (1+R^2)^(1/2)
// note: at the oriign we apply Evan's method for the terms wich present 1/r
// note: for r=s, we use a different RHS sicne some of the terms would diverge. Computing a limit "by hand", we can put = 0 the problematic terms

double wave_eq_compactified_PI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    double B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3*s*s;
    double H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    double O = pow(r,4)-pow(s,4);
    
    return (2*M*r*s*s*fields_vect[1][ind_space]/A/B
            +3*O*O*pow(s,4)*fields_vect[2][ind_space]/(pow(A,2)*B)
            +D*Dx[0](fields_vect[0],ind_space,dx)/H
            -N*r*r*s*s*Dx[0](fields_vect[1],ind_space,dx)/A/B
            -pow(s,8)/A/B* 3* (pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3))
            );
}

double wave_eq_compactified_PHI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    double B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3*s*s;
    double H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    double O = pow(r,4)-pow(s,4);
    
    return ( 2*C*F*G*r*(r+s)*fields_vect[1][ind_space] /(A*B)
            +(3*F*G*L*L*fields_vect[2][ind_space] )/(A*A*B)
            +F*s*s*Dx[0](fields_vect[0],ind_space,dx)/H
            +D*Dx[0](fields_vect[1],ind_space,dx)/H
            );
}

double wave_eq_compactified_phi_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    return(fields_vect[0][ind_space]);
}
// ----------- // MODEL 1 // ----------- //


double model1_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +param[0]*(pow(fields_vect[1][ind_space],2)-pow(fields_vect[0][ind_space],2)) );
}


// ----------- // MODEL 3 // ----------- //


double model3_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[2][ind_space] + param[0] * fields_vect[5][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2))
            );
}

double model3_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[5][ind_space] - param[0] * fields_vect[2][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2))
            );
}

double model3_PHI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[0],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double model3_PHI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[3],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double model3_phi1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[0][ind_space]);
}

double model3_phi2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[3][ind_space]);
}

// ---------- INITIAL DATA AND INITIALIZATION OF VECTORS --------------- //

double initial_line(double x,double init_param) 
{
    return(x);
}

double initial_unity(double x,double init_param) 
{
    return(1.);
}

double initial_null(double x,double init_param) 
{
    double a = 0.;
    return(a);
}


double initial_gauss_PI(double x,double init_param)
{
    return( -init_param*exp(-pow((x*2/3),2)) );
}

double initial_gauss_PHI_compactified(double x,double init_param)
{
    double s = 1;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    if(x!=s)
    {
        return( 2*init_param *pow(R,2)*exp(-pow(R,2))-init_param *exp(-pow(R,2)));
    }
    else
    {
    return(0);
    }
}

double initial_gauss_phi_compactified(double x,double init_param)
{
    double s = 1;
    double dev_std = 4;
    double rate_of_square = pow(dev_std*x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    if(x!=s)
    {
        return( -init_param *R*exp(-pow(dev_std*R,2)));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PHI(double x,double init_param)
{
    return(+2*x*init_param*exp(-pow(x,2)));
}

// ------ hyperboloidal compactification and Chi function rescaling  ------ //

double initial_gauss_PI_compactified_Chi(double x,double init_param)
{
    double s = 2.5;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    if(x!=s)
    {
        return( -init_param * Chi*exp(-pow(R,2)));
    }
    else
    {
    return(0);
    }
}
double initial_gauss_PHI_compactified_Chi(double x,double init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;

    if(x!=s)
    {
        return( ChiPrime*init_param *exp(-pow(dev_std*R,2))+Chi*(init_param*(-2*dev_std*dev_std*R)*exp(-pow(dev_std*R,2))) );
    }
    else
    {
    return(0);
    }
}

double initial_gauss_PHI_compactified_Chi1(double x,double init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;

    if(x!=s)
    {
        return( init_param*(-2*dev_std*dev_std*R)*exp(-pow(dev_std*R,2)) );
    }
    else
    {
    return(0);
    }
}

double initial_gauss_phi_compactified_Chi(double x,double init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);

    if(x!=s)
    {
        return( init_param*Chi*exp(-pow(dev_std*R,2)));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_phi_compactified_Chi1(double x,double init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);

    if(x!=s)
    {
        return( init_param*exp(-pow(dev_std*R,2)));
    }
    else
    {
        return(0);
    }
}


// ------  ------  ------  ------  ------  ------  ------  ------  ------  ------

double initial_gauss_function(double x,double init_param) 
{
    
    double a = exp(-pow(x,2));
    return(a);
}

double initial_parab_function(double x,double init_param)
{
    double a = -2*pow(x,2);
    return(a);
}

double line(double x,double t)
{
    return(4*t);
}


// initialize the fields at time zero
vector<vector<double>> initialize_fields(double d_min,double d_max,double dx,std::vector<double(*)(double,double)> funcs,double param_ic,int gl, int gr,int ord)
{
    
     if (ord>gl)
    {
        gl = ord;
    }
    if(ord>gr)
    {
        gr = ord;
    }
    int N = funcs.size();
    int S =  int( ((d_max+dx*double(gr))-(d_min-dx*double(gl)) +dx/2 ) / dx) + 1;
    //cout<<S<<endl;
    std::vector<std::vector<double>> new_vect(N, vector<double>(S));
    
    
    for(int n=0;n<N;n++)
    {
        for(double j=gl;j < S-gr ;j=j+1)
        {
            double val = d_min+dx*(j-gl);
            //cout<< j <<"\n";
            new_vect[n][j] = funcs[n](val,param_ic);
        }
    }
    return(new_vect);
}


void init_func(std::vector<double> &func_vect,double d_min,double d_max,double dx,double(*func)(double,double),double t)
{
    for(double j=d_min;j<d_max+dx/2.;j+=dx)
    {
        func_vect.push_back(func(j,t));
    }
}

// --------- BOUNDARY CONDITIONS --------- //

// adv eq 
void adv_boundaries_TEM(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //same of the bulk
    fields_vect_new[j].insert(fields_vect_new[j].begin(),((-4*fields_vect_old[0][gl]+7.*fields_vect_old[0][gl+1]-4.*fields_vect_old[0][gl+2]+fields_vect_old[0][gl+3])/2/dx+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt)));
    //boundary conditions
    fields_vect_new[j].push_back((-(4*fields_vect_old[0][last_ind]-7.*fields_vect_old[0][last_ind-1]+4.*fields_vect_old[0][last_ind-2]-fields_vect_old[0][last_ind-3]) +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
    //fields_vect_new[j].push_back(((5.*fields_vect_old[j][last_ind]-11.*fields_vect_old[j][last_ind-1]+10.*fields_vect_old[j][last_ind-2]-5*    fields_vect_old[j][last_ind-3]+fields_vect_old[j][last_ind-4])/2/dx +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
}

void adv_boundaries_right(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //same of the bulk
    fields_vect_new[j][gl] = Dx[0](fields_vect_old[0],gl,dx) + artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt);
    //boundary conditions
    fields_vect_new[j][last_ind] = (-param[0]*Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

void adv_boundaries_left(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx)+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    //fields_vect_new[j].insert(fields_vect_new[j].begin(),(-3./2.*fields_vect_old[0][gl]+2.*fields_vect_old[0][gl+1]-1*2.*fields_vect_old[0][gl+2])+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    
    

    // same of the bulk
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
    //fields_vect_new[j].push_back((-(3./2.*fields_vect_old[0][last_ind]-2.*fields_vect_old[0][last_ind-1]+1/2.*fields_vect_old[0][last_ind-2] )+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
}



void adv_boundaries_periodic(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //left boundary
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl+1]-fields_vect_old[0][last_ind])/2/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind-1])/2/dx +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
    
}

void adv_boundaries_periodic_backward_der(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    //left boundary
    
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind])/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));   
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][last_ind]-fields_vect_old[0][last_ind-1])/dx +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));

}


// wave equation

void refl_abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( (Dx[0](fields_vect_old[0],gl,dx)-Dx[0](fields_vect_old[1],gl,dx))/2 +artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt)));
    fields_vect_new[j].push_back((-(Dx[0](fields_vect_old[0],last_ind,dx)-Dx[0](fields_vect_old[1],last_ind,dx))/2+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
}



void radiative_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    int ind_space = gl;
    // left boundary conditions
    fields_vect_new[j][gl] = ( Dx[0](fields_vect_old[0],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));

    //fields_vect_new[j].insert(fields_vect_new[j].begin(), Dx[0](fields_vect_old[0],gl,dx));
    //fields_vect_new[j].insert(fields_vect_new[j].begin(),( (-3/2*fields_vect_old[0][gl]+2*fields_vect_old[0][gl+1]-1/2*fields_vect_old[0][gl+2])/dx +artificial_diss(epsilon,ord,fields_vect_old,0,gl,dx,dt)));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));

    //fields_vect_new[j].push_back(-Dx[0](fields_vect_old[0],last_ind,dx));
    //fields_vect_new[j].push_back(- ((3/2*fields_vect_old[0][last_ind]-2*fields_vect_old[0][last_ind-1]+1/2*fields_vect_old[0][last_ind-2])/dx +artificial_diss(epsilon,ord,fields_vect_old,0,last_ind,dx,dt)));
}

void radiative_outer_boundaries_PI_we(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = (3.* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) );
    fields_vect_new[0][ind_space] = ( derivative+artificial_diss(epsilon,ord,fields_vect_old,0,gl,dx,dt));

    // right boundary conditions
    fields_vect_new[0][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,0,last_ind,dx,dt));

    //fields_vect_new[j].push_back(-Dx[0](fields_vect_old[0],last_ind,dx));
    //fields_vect_new[j].push_back(- ((3/2*fields_vect_old[0][last_ind]-2*fields_vect_old[0][last_ind-1]+1/2*fields_vect_old[0][last_ind-2])/dx +artificial_diss(epsilon,ord,fields_vect_old,0,last_ind,dx,dt)));
}


void abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( Dx[0](fields_vect_old[0],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt)));
    fields_vect_new[j].push_back((-Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
}

void abs_boundaries_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j].insert(fields_vect_new[j].begin(), (Dx[0](fields_vect_old[1],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt)));
    // right boundary conditions
    fields_vect_new[j].push_back(( -Dx[0](fields_vect_old[1],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt)));
}

   
void no_boundary_conditions_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

void no_boundary_conditions_phi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[0][gl]+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[0][last_ind]+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

// compactified and rescaled wave equation//

// R rescaling
void no_boundary_conditions_PI_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    
        
    double r = dmin;
    double s = param[0];
    
    double rate_of_square = pow(r,2)/pow(s,2);
    
    //w we introduce a variable for 1/(R'(1-H'^2))
    double coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variavle for H'/(R'*(1-H'^2))
    double coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
    
     
    double left_value = (-coefficient2*Dx[0](fields_vect_old[0],ind_space_left,dx)
                        -coefficient1*Dx[0](fields_vect_old[1],ind_space_left,dx));
    
    
    r = dmax;
    s = param[0];
    
    rate_of_square = pow(r,2)/pow(s,2);
    //w we introduce a variable for 1/(R'(1-H'^2))
    coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variavle for H'/(R'*(1-H'^2))
    coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
    
    
    double right_value = (-coefficient2*Dx[0](fields_vect_old[0],last_ind,dx)
                        -coefficient1*Dx[0](fields_vect_old[1],last_ind,dx));
     
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));    
}

// R rescaling
void no_boundary_conditions_PHI_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    
    double r = dmin;
    double s = param[0];
    
    double rate_of_square = pow(r,2)/pow(s,2);
    
    //w we introduce a variable for 1/(R'(1-H'^2))
    double coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variavle for H'/(R'*(1-H'^2))
    double coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
    
     
    double left_value = (-coefficient2*Dx[0](fields_vect_old[1],ind_space_left,dx)
                        -coefficient1*Dx[0](fields_vect_old[0],ind_space_left,dx));
    
    
    r = dmax;
    
    rate_of_square = pow(r,2)/pow(s,2);
    //w we introduce a variable for 1/(R'(1-H'^2))
    coefficient1 = (1+rate_of_square) / (1+4*rate_of_square-pow(rate_of_square,2));
    
    // we introduce a variavle for H'/(R'*(1-H'^2))
    coefficient2 = rate_of_square*(3-rate_of_square)/(1+4*rate_of_square-pow(rate_of_square,2));
    
    
    double right_value = (-coefficient2*Dx[0](fields_vect_old[1],last_ind,dx)
                        -coefficient1*Dx[0](fields_vect_old[0],last_ind,dx));
     
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));    
}

// R rescaling
void no_boundary_conditions_phi_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    
    double left_value = -fields_vect_old[0][ind_space_left];
    double right_value = -fields_vect_old[0][last_ind];
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt);    
}


// CHI RESCALING //


void no_boundary_conditions_PI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
        
    double r = dmin;
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    double B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3*s*s;
    double H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    double O = pow(r,4)-pow(s,4);
    
    
    double left_value = 2*M*r*s*s*fields_vect_old[1][ind_space_left]/A/B
            +3*O*O*pow(s,4)*fields_vect_old[2][ind_space_left]/(pow(A,2)*B)
            +D*Dx[0](fields_vect_old[0],ind_space_left,dx)/H
            -N*r*r*s*s*Dx[0](fields_vect_old[1],ind_space_left,dx)/A/B
            -pow(s,8)/A/B* 3* (pow((r+dx),2)*fields_vect_old[1][ind_space_left+1] - pow((r-dx),2)*fields_vect_old[1][ind_space_left-1] )/(pow((r+dx),3)-pow((r-dx),3));
    
    r = dmax;
    
    A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    C = r- s;
    D = pow(r,4) - 3*r*r*s*s;
    F = pow(r,2) + s*s;
    G = pow(r,2)-3*s*s;
    H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    L = pow(r,3)*s - r* pow(s,3);
    M = pow(r,4) + r*r * s*s - pow(s,4);
    N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    O = pow(r,4)-pow(s,4);
    
    double right_value =2*M*r*s*s*fields_vect_old[1][ind_space_left]/A/B
            +3*O*O*pow(s,4)*fields_vect_old[2][last_ind]/(pow(A,2)*B)
            +D*Dx[0](fields_vect_old[0],last_ind,dx)/H
            -N*r*r*s*s*Dx[0](fields_vect_old[1],last_ind,dx)/A/B
            -pow(s,8)/A/B* 3* (pow((r+dx),2)*fields_vect_old[1][last_ind+1] - pow((r-dx),2)*fields_vect_old[1][last_ind-1] )/(pow((r+dx),3)-pow((r-dx),3));

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));   
}

void no_boundary_conditions_PHI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
        
    double r = dmin;
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    double B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3*s*s;
    double H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    double O = pow(r,4)-pow(s,4);
    //cout<<" A B C D F G L"<<A<<";"<<B<<";"<<C<<";"<<D<<";"<<F<<";"<<G<<";"<<L<<endl;

    
    double left_value =  2*C*F*G*r*(r+s)*fields_vect_old[1][ind_space_left] /(A*B)
            +(3*F*G*L*L*fields_vect_old[2][ind_space_left] )/(A*A*B)
            +F*s*s*Dx[0](fields_vect_old[0],ind_space_left,dx)/H
            +D*Dx[0](fields_vect_old[1],ind_space_left,dx)/H;
    
    r = dmax;
    
    A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    B = pow(r,4) - 4 * pow(r,2)* s*s - pow(s,4);
    C = r- s;
    D = pow(r,4) - 3*r*r*s*s;
    F = pow(r,2) + s*s;
    G = pow(r,2)-3*s*s;
    H = -pow(r,4) + 4*pow(r,2)*s*s + pow(s,4);
    L = pow(r,3)*s - r* pow(s,3);
    M = pow(r,4) + r*r * s*s - pow(s,4);
    N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1 + s*s);
    O = pow(r,4)-pow(s,4);
    
    double right_value =  2*C*F*G*r*(r+s)*fields_vect_old[1][last_ind] /(A*B)
            +(3*F*G*L*L*fields_vect_old[2][last_ind] )/(A*A*B)
            +F*s*s*Dx[0](fields_vect_old[0],last_ind,dx)/H
            +D*Dx[0](fields_vect_old[1],last_ind,dx)/H;
     
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));   
}

void no_boundary_conditions_phi_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
        
    double r = dmin;   
    
    double left_value = fields_vect_old[0][ind_space_left];
    
    r = dmax;

    double right_value = fields_vect_old[0][last_ind];
     
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value+artificial_diss(epsilon,ord,fields_vect_old,j,ind_space_left,dx,dt));
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));   
}
// ------------- // MODEL 1 // ------------- //


void radiative_outer_boundaries_PI_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    
    
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3.* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +(pow(fields_vect_old[1][ind_space],2)-pow(fields_vect_old[0][ind_space],2) );
    
    fields_vect_new[j][ind_space] = ( derivative+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));

    
    
    // right boundary conditions
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));

    
}



// ------------- // MODEL 3 // ------------- //


void radiative_outer_boundaries_PI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))      
                        + (fields_vect_old[2][ind_space] + param[0] * fields_vect_old[5][ind_space])/pow(param[0],2)*
                        (pow(fields_vect_old[1][ind_space],2)+pow(fields_vect_old[4][ind_space],2)-pow(fields_vect_old[0][ind_space],2)-pow(fields_vect_old[3][ind_space],2));
    
    fields_vect_new[j][ind_space] = ( derivative+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)-fields_vect_old[0][last_ind]/x +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

void radiative_outer_boundaries_PI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3* ( pow((x+dx),2)*fields_vect_old[4][ind_space+1] - pow((x-dx),2)*fields_vect_old[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) 
                        + (fields_vect_old[5][ind_space] - param[0] * fields_vect_old[2][ind_space])/pow(param[0],2)*
                        (pow(fields_vect_old[1][ind_space],2)+pow(fields_vect_old[4][ind_space],2)-pow(fields_vect_old[0][ind_space],2)- 
                        pow(fields_vect_old[3][ind_space],2));
    
                        
                        
    fields_vect_new[j][ind_space] = ( derivative+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[3],last_ind,dx)-fields_vect_old[3][last_ind]/x +artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}


void no_boundary_conditions_PHI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[0],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

void no_boundary_conditions_PHI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[3],gl,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[3],last_ind,dx)+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}


void no_boundary_conditions_phi1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[0][gl]+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[0][last_ind]+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

void no_boundary_conditions_phi2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[3][gl]+artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt));
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[3][last_ind]+artificial_diss(epsilon,ord,fields_vect_old,j,last_ind,dx,dt));
}

//-------------- DIFFERENTIAL OPERATORS -------------//
    
// simple function to get the derivative of a function, f=g'
// It uses boundary condition!
// input:primitive function (g),discretization step


    
// fourth order central finite difference expression for a first spatial derivative
double first_der_fourth_order_centered(std::vector<double> vector,int i,double dx)
{
    return((vector[i-2]/12.-2./3.*vector[i-1]+2./3.*vector[i+1]-vector[i+2]/12.)/dx);
}
    
// fourth order central finite difference expression for a first spatial derivative
double first_der_second_order_centered(std::vector<double> vector,int i,double dx)
{
    return((-vector[i-1]+vector[i+1])/2./dx);
}

double first_der_first_order_backward(std::vector<double> vector,int i,double dx)
{
    return((-vector[i-1]+vector[i])/dx);
}

double first_der_second_order_backward(std::vector<double> vector,int i,double dx)
{
    return((1/2.*vector[i-2]-2.*vector[i-1]+3./2.*vector[i])/dx);
}

double first_der_second_order_forward(std::vector<double> vector,int i,double dx)
{
    return((1/2.*vector[i-2]+2.*vector[i-1]-3./2.*vector[i])/dx);
}
    
//-------------- GHOST POINTS --------------//

void ghost_point_extrapolation_1_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    
    // attachment of the ghost points to the boundary of the function
    for (int i=gl-1; i>=0;i--)
    {
        field_vect[j][i] = (2*field_vect[j][i+1]-field_vect[j][i+2]);
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = (2*field_vect[j][i-1]-field_vect[j][i-2]);
    }
}

void ghost_point_extrapolation1(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    // computation of the ghost points via interpolation
    // we find the x vector and the relative y vector values of the subsections starting from the boundaries
    int number_of_point_int = 30;
    std::vector<double> left_side_y;
    std::vector<double> left_side_x;
    for (int i = 0;i<number_of_point_int;i++)
    {
        left_side_x.push_back(dmin+i*dx);
        left_side_y.push_back(field_vect[j][i+gl]);
    }
        

    std::vector<double> right_side_x;
    std::vector<double> right_side_y;
    for (int i = number_of_point_int-1;i>=0;i--)
    {
        right_side_x.push_back(dmax-i*dx);
        //cout<<"x:"<<(dmax-i*dx);
        right_side_y.push_back(field_vect[j][field_vect[j].size()-1-i-gr]);
        //cout<<" y:"<<(field_vect[j][field_vect[j].size()-1-i-gr]);

    }
        
    // interpolation 
    tk::spline left(left_side_x,left_side_y);
    tk::spline right(right_side_x,right_side_y);
        
    // attachment of the ghost points to the boundary of the function
        
    for (int i=gl-1; i>=0;i--)
    {
        field_vect[j][i] = left(dmin-dx*i-dx);
    }
        
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
    int counter = 1;
    field_vect[j][i] = right(dmax+dx*counter);
    counter = counter +1;
    }
}


void ghost_point_extrapolation2_TEM(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int last_ind = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
        
    for (int i=gl-1; i>=0;i--)
    {
        field_vect[j][i] = (field_vect[j].begin(),(field_vect[j][0]-(-4*field_vect[j][0]+7.*field_vect[j][1]-4.*field_vect[j][2]+field_vect[j][3])/2));
    }
        
    for (int i=0; i<gr;i++)
    {
        field_vect[j].push_back(field_vect[j][last_ind]+(4*field_vect[j][last_ind]-7.*field_vect[j][last_ind-1]+4.*field_vect[j][last_ind-2]-field_vect[j][last_ind-3])/2);
    }
}

void ghost_point_extrapolation_2_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int last_ind = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    for (int i=gl-1; i>=0;i--)
    {
        field_vect[j][i] = field_vect[j][i+1]-(-3./2.*field_vect[j][i+1]+2.*field_vect[j][i+2]-1./2.*field_vect[j][i+3])+1./2.*(2.*field_vect[j][i+1]-5.*field_vect[j][i+2]+4.*field_vect[j][i+3]-field_vect[j][i+3]);
    }
        
    for (int i=last_ind-gr; i<last_ind;i++)
    {
        cout<<i<<endl;
        field_vect[j][i] = field_vect[j][i-1]+(3./2.*field_vect[j][i-1]-2.*field_vect[j][i-2]+1./2.*field_vect[j][i-3])+1./2.*(2.*field_vect[j][i-1]-5.*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4]);
    }
}

void ghost_point_extrapolation_4_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    ;
    // attachment of the ghost points to the boundary of the function
    for (int i=gl-1; i>=0;i--)
    {
        
        field_vect[j][i] = field_vect[j][i+1]-(-25./12.*field_vect[j][i+1]+4.*field_vect[j][i+2]-3.*field_vect[j][i+3]+4./3.*field_vect[j][i+4]-1./4.*field_vect[j][i+5])+1./2.*(2.*field_vect[j][i+1]-5.*field_vect[j][i+2]+4.*field_vect[j][i+3]-field_vect[j][i+4]);
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = field_vect[j][i-1]+(+25./12.*field_vect[j][i-1]-4.*field_vect[j][i-2]+3.*field_vect[j][i-3]-4./3.*field_vect[j][i-4]+1./4.*field_vect[j][i-5]+1./2.*(2.*field_vect[j][i-1]-5.*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4]));
    }
}

void ghost_point_extrapolation_6_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    ;
    // attachment of the ghost points to the boundary of the function
    for (int i=0; i<gl;i++)
    {
        field_vect[j].insert(field_vect[j].begin(),(field_vect[j][0]-(-49./20.*field_vect[j][0]+6.*field_vect[j][1]-15./2.*field_vect[j][2]+20./3.*field_vect[j][3]-15./4.*field_vect[j][4]+6./5.*field_vect[j][5]-1./6.*field_vect[j][6])));
    }
        
    for (int i=0; i<gr;i++)
    {
        int last_ind = field_vect[j].size()-1;
        field_vect[j].push_back(field_vect[j][last_ind]+(49./20.*field_vect[j][last_ind]-6.*field_vect[j][last_ind-1]+15./2.*field_vect[j][last_ind-2]-20./3.*field_vect[j][last_ind-3]+15./4.*field_vect[j][last_ind-4]-6./5.*field_vect[j][last_ind-5]+1./6.*field_vect[j][last_ind-6]));
    }
}

void ghost_point_extrapolation3(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
        
    // attachment of the ghost points to the boundary of the function
        
    for (int i=0; i<gl;i++)
    {
        field_vect[j].insert(field_vect[j].begin(),field_vect[j][0]);
    }
        
    for (int i=field_vect.size()-1; i>field_vect.size()-1-gl;i--)
    {
    field_vect[j].push_back(field_vect[j][field_vect.size()-1]);
    }
}

void ghost_point_extrapolation_1_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = field_vect[j][i] = (2*field_vect[j][i-1]-field_vect[j][i-2]);
    }
    
}

void ghost_point_extrapolation_4_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = +1./24.*(109*field_vect[j][i-1]-194.*field_vect[j][i-2]+170.*field_vect[j][i-3]-76.*field_vect[j][i-4]+17.*field_vect[j][i-5]-2.*field_vect[j][i-6]);
    }
}

void ghost_point_extrapolation_1_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = field_vect[j][i-1]+(4*field_vect[j][i-1]-7.*field_vect[j][i-1-1]+4.*field_vect[j][i-1-2]-field_vect[j][i-1-3])/2;
    }
    
}

void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = 4.*field_vect[j][i-1]-6*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4];
    }
    
}

void ghost_point_extrapolation_2_ord_spline_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = field_vect[j][i-1]+(4*field_vect[j][i-1]-7.*field_vect[j][i-1-1]+4.*field_vect[j][i-1-2]-field_vect[j][i-1-3])/2;
    }
    
}

void ghost_point_extrapolation_4_ord_spherical_symmetry_rescaled(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PI and phi are odd
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = -1.*field_vect[j][2*gl-i];
        }
    }
    
    for (int i=S-gr; i<S;i++)
    {
        field_vect[j][i] = (field_vect[j][i-1]+(+25./12.*field_vect[j][i-1]-4.*field_vect[j][i-2]+3.*field_vect[j][i-3]-4./3.*field_vect[j][i-4]+1./4.*field_vect[j][i-5]));
    }
}

void ghost_point_extrapolation_6_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax)
{
    ;
    // attachment of the ghost points to the boundary of the function
    if (j==0 || j==2)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j].insert(field_vect[j].begin(),field_vect[j][1+i]);
        }
    }
    
    // PHI is odd
    if (j==1)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j].insert(field_vect[j].begin(),-field_vect[j][1+i]);
        }
    }
        
    for (int i=0; i<gr;i++)
    {
        int last_ind = field_vect[j].size()-1;
        field_vect[j].push_back(field_vect[j][last_ind]+(49./20.*field_vect[j][last_ind]-6.*field_vect[j][last_ind-1]+15./2.*field_vect[j][last_ind-2]-20./3.*field_vect[j][last_ind-3]+15./4.*field_vect[j][last_ind-4]-6./5.*field_vect[j][last_ind-5]+1./6.*field_vect[j][last_ind-6]));
    }
}


// ------------- ARTIFICIAL DISSIPATION -------------- //
double artificial_dissipation_2(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*(dt/dx)*pow(-1,ord)*(copy_fields_vect[j][i+2]-4*copy_fields_vect[j][i+1]+6*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_Husa(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon/pow(dx,1)*pow(-1,ord)/16*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
}

// to control!
double artificial_dissipation_4(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*pow(dx,3)*pow(-1,ord)/4*(copy_fields_vect[j][i+2]-4*copy_fields_vect[j][i+1]+6*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_ani(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*pow(dx,2*ord-1)*pow(-1,ord)/4.*(copy_fields_vect[j][i+2]-4*copy_fields_vect[j][i+1]+6*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4*copy_fields_vect[j][i-1]) );
}
    
// ------------- NORMS AND CONVERGENCY -------------- //

double norm(std::vector<double> &vector, double h)
{
    double norm = 0;
    for (int i=0; i<vector.size(); i++)
    {
        norm = norm + pow(vector[i],2);
    }
    norm = norm * h;
    norm = pow(norm,0.5);
    return(norm);
}

double norm_of_diff(std::vector<double> &num1, std::vector<double> &num2,double h, double(*norm) (std::vector<double> & , double ))
{
    std::vector<double> diff;
    for (int i=0; i<num2.size(); i++)
    {
        diff.push_back(num1[i]-num2[i]);
    }
    return( norm(diff,h) );
}

double conv_test(std::vector<double> &num1,std::vector<double> &num2,double(*theo_sol)(double,double),void(*init_func)(std::vector<double> &,double ,double ,double ,double(*)(double,double),double),double h1,double h2, double(*norm) (std::vector<double> & , double ),double (*norm_of_diff)(std::vector<double> &, std::vector<double> &,double , double(*) (std::vector<double> & , double )),double dmin, double dmax,double dx,double t)
{
    cout<<"--- Perfoming convergence test ---\nrate of the deltas="<<h1/h2<<endl;
    std::vector<double> theo1;
    init_func(theo1,dmin,dmax,h1,theo_sol,t);
    //cout<<"1 "<<theo1.size()<<" "<<num1.size()<<endl;
    double norm_h1 = norm_of_diff(num1,theo1,h1,norm);
    //cout<<"norm1 "<<norm_h1<<"\n";
    
    
    std::vector<double> theo2;
    init_func(theo2,dmin,dmax,h2,theo_sol,t);
    //cout<<"2 "<<theo2.size()<<" "<<num2.size()<<endl;
    double norm_h2 = norm_of_diff(num2,theo2,h2,norm);
    //cout<<"norm2 "<<norm_h2<<"\n";
    
    double base = h1/h2;
    cout<<"conv test"<<log(norm_h1/norm_h2)/log(base) <<"\n";
    return(log(norm_h1/norm_h2)/log(base)  );
}


double self_conv_test(std::vector<double> &num1,std::vector<double> &num2,std::vector<double> &num3, double h1, double h2, double(*norm) (std::vector<double> & , double ),double (*norm_of_diff)(std::vector<double> &, std::vector<double> &,double, double(*) (std::vector<double> & , double )) )
{
    cout<<"--- Perfoming self convergence test ---\n";
    std::vector<double> support12;
    for (int i=0; i<num2.size(); i = i+2)
    {
        support12.push_back(num2[i]);
    }
    //cout<<"num1 "<<num1.size() << "\nnum2 "<<num2.size()<<"\nsupport "<<support12.size()<<"\n";
    double norm_h1 = norm_of_diff(num1,support12,h1,norm);
    
    std::vector<double> support23;
    for (int i=0; i<num3.size(); i = i+2)
    {
        support23.push_back(num3[i]);
    }
    //cout<<"num3 "<<num3.size() << "\nnum2 "<<num2.size()<<"\nsupport23 "<<support23.size()<<"\n";
    double norm_h2 = norm_of_diff(num2,support23,h2,norm);
    cout<<"self conv test "<<log(norm_h1/norm_h2)/log(2)<<"\n";
    return(log(norm_h1/norm_h2)/log(2));
}




void diff_vector(std::vector<double> &diff,std::vector<double> &vec1, std::vector<double> &vec2)
{
    for(int i=0; i<vec1.size();i++)
    {
        diff.push_back(vec1[i]-vec2[i]);
    }
}

// ------------- PRINTING AND READING FUNCTIONS --------------- //
    
// Function to print on a file the data



void print_f(std::vector< std::vector<double> > fields_vect, double dmin, double dx, string name_file,string name_folder, int gl, int gr)
{
    
    ofstream myfile_print_f;
    myfile_print_f.open (name_file,ios::app);
    // headers of the columns
    myfile_print_f<<"x,";
    for (int i=0; i<fields_vect.size()-1; i++)
    {
        myfile_print_f << "field"<<to_string(i)<<",";
    }    
    myfile_print_f << "field"<<to_string(fields_vect.size()-1)<<"\n";
    
    // for every spatial point
    for (int j=0; j<fields_vect[0].size();j= j+1)
    {
        myfile_print_f << dmin+dx*double(j-gl)<<",";
        // for every fields add the relative value
        for (int i=0; i<fields_vect.size()-1; i++)
        {
            myfile_print_f << fields_vect[i][j]<<",";
        }
        myfile_print_f << fields_vect[fields_vect.size()-1][j];
        myfile_print_f<<"\n";
    
    }
    //cout<<" field0"<<"= "<< fields_vect[0][fields_vect[0].size()-1]<<endl;
    
    myfile_print_f.close();
}


void read_parameters(string name_parameters_file, double &dmin, double &dmax, double &h1, double &integration_interval,int &step_to_save,int &gl,int &gr,int &ord, vector<double> &epsilon,vector<double> &parameters_ic_vector, vector<double> &parameters)
{
    // we read the input parameter from an external file
    ifstream input_file;
    input_file.open(name_parameters_file);
    
    input_file.ignore(256,' ');
    input_file >> dmin;
    
    input_file.ignore(256,' ');
    input_file >> dmax;
    
    input_file.ignore(256,' ');
    input_file >> h1;
    
    input_file.ignore(256,' ');
    input_file >> integration_interval;
    
    input_file.ignore(256,' ');
    input_file >> step_to_save;
    
    input_file.ignore(256,' ');
    input_file >> gl;
    
    input_file.ignore(256,' ');
    input_file >> gr;
    
    input_file.ignore(256,' ');
    input_file >> ord;
    
    
    // epsilon vector
    std::string e;
    input_file.ignore(256,' ');
    getline (input_file,e);
    
    std::stringstream iss (e);
    double number;
    while ( iss >> number )
    {
        epsilon.push_back( number );
    }
    
    // initial parameters vector
    std::string ic_par;
    input_file.ignore(256,' ');
    getline (input_file,ic_par);
    std::stringstream iss1 (ic_par);
    while ( iss1 >> number )
    {
        parameters_ic_vector.push_back( number );
    }
    
    // parameters vector
    std::string par;
    input_file.ignore(256,' ');
    getline (input_file,par);
    std::stringstream iss2 (par);
    while ( iss2 >> number )
    {
        parameters.push_back( number );
    }
    
    input_file.close();
}
