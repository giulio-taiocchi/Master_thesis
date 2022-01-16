#include <math.h>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include "spline.h"
using namespace std;

typedef vector<vector<vector< vector<double> >>> grid;
typedef vector<vector<vector< vector<double> >>> fields_vector;


typedef void (*print_function)(fields_vector fields_vect, vector<double> &dmin,vector<double> &dmax, vector<double> dx, string name_file,string name_folder, int gl, int gr,MPI_Status status, int totalnodes, int mynode,MPI_Request request);


typedef double (*artificial_dissipation_function)(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,vector<double> dx,double dt);

typedef std::vector<double (*)(std::vector<double>,int,double)>  derivative_vector;

typedef void(*ghost_point_extrapolation_function)(fields_vector &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax);

typedef double (*evolution_function)(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl);

typedef void(*boundary_conditions_function)(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo);

typedef void(*communication_function)(std::vector< std::vector<double> > &fields_vect,int j,int index_dmax_local,int index_dmin_local,int nitems,int mynode,int totalnodes,MPI_Status status);

typedef vector<vector<double>>(*one_step_function)(vector<vector<vector<double>>> fields_vect,vector<double> &dmin,vector<double> &dmax,vector<double> dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord,MPI_Status status, int totalnodes, int mynode, communication_function communication)   ;

typedef void(*method_of_line_function)(vector<vector<vector<double>>> fields_vect,one_step_function one_step,int dim, vector<double> dx, std::vector<double> param, double dt, double interval,vector<double> dmin,vector<double> dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path ,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication);

typedef fields_vector(*initialization_fields_function)(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> &param_ic,int gl, int gr,int ord,int dim, grid Grid);

typedef grid(*initialization_grid_function)(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim);



void multiple_parameters_run(std::vector<double>& parameters_ic_vector, std::vector<double(*)(double, double)>& initial_conditions, initialization_grid_function initialize_grid, initialization_fields_function initialize_fields,int dim, vector<double> &dmin, vector<double> &dmax, double h1, double h2, double h3, double dt1, double dt2, double dt3, double integration_interval,int step_to_save, derivative_vector & Dx ,std::vector< evolution_function > & R_vector, std::vector< boundary_conditions_function > & b_func,  std::vector<double> & parameters, one_step_function one_step, int gl,int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2, vector<double> epsilon1,print_function print_f, string file_path, method_of_line_function MOL_RK4,int ord,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication)
{
    int dim = 3;
    ofstream myfile2;
    //myfile2.open (file_path+"name_of_file");
    //myfile2<<"names\n";
    //myfile2.close();
    
    for(int e=0;e<epsilon1.size();e++)
    {
        for (int l=0;l<parameters_ic_vector.size();l++)
        {
            //cout<<"processor :"<<mynode<<"\nepsilon = "<<epsilon1[e]<<" amplitude = "<<parameters_ic_vector[l]<<endl;
            
            // initialize fields
            
            
            grid Grid = = initialize_grid(dmin,dmax,dx, gl,gr,ord,dim);
            fields_vector fields_vect1 = initialize_fields(dmin,dmax,h1,initial_conditions,parameters_ic_vector[l],gl,gr,ord);     
            fields_vector fields_vect2 = initialize_fields(dmin,dmax,h2,initial_conditions,parameters_ic_vector[l],gl,gr,ord);
            fields_vector fields_vect3 = initialize_fields(dmin,dmax,h3,initial_conditions,parameters_ic_vector[l],gl,gr,ord);
            // setting the output variables
            //cout<<"\n lun vec \n"<<endl;
            string name_file = "processor_"+to_string(mynode)+"_ampl_"+to_string(parameters_ic_vector[l])+"_eps"+to_string(epsilon1[e])+"_dx_"+to_string(h1[0])+"steps"+to_string(step_to_save)+"last_time"+to_string(integration_interval)+".csv";
            
            string complete_path = file_path+name_file;
            ofstream myfile;
            
            myfile.open (complete_path);
            myfile.close();
            
            MOL_RK4(fields_vect1,one_step,h1,parameters,dt1, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path,status,totalnodes,mynode, request,communication);
            MOL_RK4(fields_vect2,one_step,h2,parameters,dt2, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path,status,totalnodes,mynode, request, communication);
            MOL_RK4(fields_vect3,one_step,h3,parameters,dt3, integration_interval,dmin,dmax,R_vector,b_func,step_to_save,print_f,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon1[e],ord,Dx,complete_path,status,totalnodes,mynode, request, communication);
        }
    }
}
    

void MOL_RK4(fields_vector fields_vect,one_step_function one_step,int dim, vector<double> dx, std::vector<double> param, double dt, double interval,    vector<double> dmin,    vector<double> dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication)
{   
    
    cout<<"\nProcessor :"<<mynode<<endl<<"--- Method of lines called ---\ndx = "<<dx[0]<<"\ndt = "<<dt<<"\nDomain = ["<<dmin[0]<<","<<dmax[0]<<"]\nlast time objective :"<<interval<<"\n";
    
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
    //extrapolation of the ghost point for the initial data
    int N = fields_vect.size();
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, 0,dx,dt,j,gl,gr,dmin,dmax);
    }
    for (double t=0;t<interval+dt/1.5;t=double(t)+double(dt)) 
    {        
        
        // we print to file only for some time
        if ((counter%m)==0)
        {
            
            //ofstream myfile;
            //myfile.open (file_path, ios::app);
            //myfile<<t<<"\n";
            //myfile.close();
            //cout<<"print the fields at time"<<t<<endl;
            print_f(fields_vect,dmin,dmax,dx,file_path,file_path,gl,gr,status,totalnodes,mynode,request); // the print_f function is called
            //cout<<t<<endl;
        }
        fields_vect = one_step(fields_vect,dmin,dmax,dx,param,dt,R_vect,bc,t,gl,gr,ghost_point_extrapolation,artificial_diss_2,epsilon,Dx,ord
                        ,status,totalnodes,mynode,communication);
        
        counter += 1;
    last = t;
    }
    
    //cout<<"last time included:"<<last<<endl;
    
}


vector<vector<double>> onestep_RK4_1(vector<vector<vector<double>>> fields_vect,vector<double> &dmin,vector<double> &dmax,vector<double> dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord,MPI_Status status, int totalnodes, int mynode, communication_function communication)  
{
    //cout<<" processor: "<<mynode<<"time "<<t<<endl;
    int nitems=2;
    int N = fields_vect.size();
    int S =  fields_vect[0].size();
    
    vector<vector<double>> k1(N, vector<double>(S)) ;
    std::vector<std::vector<double>> support_k1(N, vector<double>(S));
    std::vector<std::vector<double>> k2(N, vector<double>(S));
    std::vector<std::vector<double>> support_k2(N, vector<double>(S));
    std::vector<std::vector<double>> k3(N, vector<double>(S));
    std::vector<std::vector<double>> support_k3(N, vector<double>(S));
    std::vector<std::vector<double>> k4(N, vector<double>(S));
    std::vector<std::vector<double>> support_k4(N, vector<double>(S));
       
    double index_dmin_local,index_dmax_local;
        
    // defining the first and last spatial index of the domain, the actual process will work inside this range
        
    index_dmin_local = (gl+1) + mynode * int((S-2-gl-gr)/totalnodes);
    //cout<<"resolution: "<<dx<<" processor: "<<mynode<<" index_dmin_local: "<<index_dmin_local<<endl;
    if (mynode==totalnodes-1)
    {
        index_dmax_local =  S-1-gr;
    }
    else
    {
        index_dmax_local = (gl+1) + (mynode+1) * int((S-2-gl-gr)/totalnodes);
    }
    //cout<<"resolution: "<<dx<<"processor: "<<mynode<<" index_dmax_local: "<<index_dmax_local<<endl;

    
    // populate the ghost zone of the fields, we need them to calculate k1
    /*
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }*/
    
    // k1 building
    for (int j=0; j <N; j++)
    {
        // evualuating the "bulk" of k1
        // we have to consider the physical domain of the fields (so exclude the GP) and then exclude the boundaries value (+1 and -1)
        // we are also dividing the spatial domain in subdomains, in which each processor run in parallel
        
        // calculating the k1 function (just in the interval associated to the actual processor)
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k1[j][i] = (evo[j](j,i,fields_vect,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl));  
                            //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        communication(k1,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
    
        
        // evaluating the boundary of k1
        bc[j](k1,fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        //k1[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt);
        //k1[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt);
        // GP extrapolation for k1 vector
        ghost_point_extrapolation(k1, t,dx,dt,j,gl,gr,dmin,dmax);

    // computing the argument for the next coefficient k2
        if(mynode==0)
        {
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];
            }
        }
        if(mynode==totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];
            }
        }
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];
            }
        }
    }
    
    
    
    // k2 building
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k2[j][i] = evo[j](j,i,support_k1,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt/2,gl);
            //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2));
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(k2,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
        //boundary conditions update
        bc[j](k2,support_k1, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        // artificial dissipation
        //k2[j][gl] += artificial_diss(epsilon,ord,support_k1,j,gl,dx,dt/2);
        //k2[j][gr] += artificial_diss(epsilon,ord,support_k1,j,gr,dx,dt/2);
        
        
        ghost_point_extrapolation(k2, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax);
        
        if(mynode==0)
        {
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];
            }  
        }
        if(mynode==totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];
            }
        }
        
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];
            }
        }
    }
    
    // k3 building
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k3[j][i] = evo[j](j,i,support_k2,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt,gl);
            //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(k3,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        bc[j](k3,support_k2, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        //k3[j][gl] += artificial_diss(epsilon,ord,support_k2,j,gl,dx,dt/2);
        //k3[j][gr] += artificial_diss(epsilon,ord,support_k2,j,gr,dx,dt/2);
        
        ghost_point_extrapolation(k3, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax);
        if(mynode==0)
        {
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];
            }
        }
        
        if(mynode==totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];
            }
        }
        
        if(mynode!=0 && mynode!=totalnodes)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];
            }
        }
        
    }
   
   
    // k4 building    
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k4[j][i] = evo[j](j,i,support_k3,dx,dmin,param,t+dt,Dx,artificial_diss,epsilon,ord,dt/2,gl);
            //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
        }
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        communication(k4,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
        bc[j](k4,support_k3, t+dt,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        //k4[j][gl] += artificial_diss(epsilon,ord,support_k3,j,gl,dx,dt);
        //k4[j][gr] += artificial_diss(epsilon,ord,support_k3,j,gr,dx,dt);
        
        ghost_point_extrapolation(k4, t+dt,dx,dt,j,gl,gr,dmin,dmax);
    }
    // we create a new vector that contains all the new fields. It is a support vector that will be swapped with the old one
    
    std::vector<std::vector<double>> new_fields_vect(N, vector<double>(S));

    for (int j=0; j <N; j++)
    {
        if(mynode==0)
        {
             for (int i=gl;i<index_dmax_local;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.
                +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
            ghost_point_extrapolation(new_fields_vect,t+dt,dx,dt,j,gl,gr,dmin,dmax);
        }
        
        if(mynode==totalnodes-1)
        {
            for (int i=index_dmin_local;i<S;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.
                                        +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
            ghost_point_extrapolation(new_fields_vect,t+dt,dx,dt,j,gl,gr,dmin,dmax);
        }
        
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local;i<index_dmax_local;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2*k2[j][i]+2*k3[j][i]+k4[j][i])/6.
                                        +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
        }
        communication(new_fields_vect,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
    }   
    

    // populate the ghost zone of the fields, we need them to calculate k1 in the next time step
    /*for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(new_fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }
    //cout<<"old "<<fields_vect.size()<<"new "<<new_fields_vect.size()<<endl;
    */
    return(new_fields_vect);
}



// ------------ TYPES OF D.E. --------------- //



// ----------- // ADVECTION EQUATION // ----------- //



// advection equation, central finite difference, 2nd order method
double advection_eq_left_going(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    {
        return (param[0]*Dx[0](fields_vect[0],ind_space,dx));
    }
}

double advection_eq_right_going(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    {
        return (param[0]*(-1.)*Dx[0](fields_vect[0],ind_space,dx));
    }
    
}

// ----------- // WAVE EQUATION // ----------- //


double wave_eq_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    return (Dx[0](fields_vect[1],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[1],ind_space,dx));

}
    
double wave_eq_PHI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[0],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double wave_eq_phi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[0][ind_space]);
}


//  wave equation in spherical symmetry
double wave_eq_spherical_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt,int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    // Dx[0](fields_vect[1],ind_space,dx) + 2./x * fields_vect[1][ind_space]
    // 3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) 
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) );
}

// R RESCALING HYPERBOLOIDAL foliation //
// here are reported the evolution equations for the auxiliary function of the rescaled (by R) and compactified wave equation
double wave_eq_compactified_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
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

double wave_eq_compactified_PHI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
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

double wave_eq_compactified_phi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    return(-fields_vect[0][ind_space]);
}


//------------------- CHI RESCALING HYPERBOLOIDAL FOLIATION -------------------//

// compactification: T = t + H(R), H'=1-1/R', R=r/(1-(r/s)^2)
// rescaling: by chi(R) = (1+R^2)^(1/2)
// note: at the oriign we apply Evan's method for the terms wich present 1/r
// note: for r=s, we use a different RHS sicne some of the terms would diverge. Computing a limit "by hand", we can put = 0 the problematic terms

double wave_eq_compactified_PI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2 + s*s);
    double B = pow(r,4) - 4. * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3.*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3*s*s;
    double H = -pow(r,4) + 4.*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1. + s*s);
    double O = pow(r,4)-pow(s,4);
    
    return (2*M*r*s*s*fields_vect[1][ind_space]/A/B
            +3*O*O*pow(s,4)*fields_vect[2][ind_space]/(pow(A,2)*B)
            +D*Dx[0](fields_vect[0],ind_space,dx)/H
            -N*r*r*s*s*Dx[0](fields_vect[1],ind_space,dx)/A/B
            -pow(s,8)/A/B* 3.* (pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3))
            );
}

double wave_eq_compactified_PHI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
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

double wave_eq_compactified_phi_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    return(fields_vect[0][ind_space]);
}
// ----------- // MODEL 1 // ----------- //


double model1_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +param[0]*(pow(fields_vect[1][ind_space],2)-pow(fields_vect[0][ind_space],2)) );
}


// ----------- // MODEL 3 // ----------- //


double model3_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[2][ind_space] + param[0] * fields_vect[5][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2))
            );
}

double model3_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[5][ind_space] - param[0] * fields_vect[2][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2))
            );
}

double model3_PHI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[0],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double model3_PHI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (Dx[0](fields_vect[3],ind_space,dx));
    //return (param.at(0)*Dx[0](fields_vect[0],ind_space,dx));
}

double model3_phi1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[0][ind_space]);
}

double model3_phi2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,vector<double> dx,vector<double> &dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //return(fields_vect[0][ind_space]);
    return(fields_vect[3][ind_space]);
}

// ---------- INITIAL DATA AND INITIALIZATION OF VECTORS --------------- //

double initial_line(vector<double> x,vector<double> init_param) 
{
    return(x);
}

double initial_unity(vector<double> x,vector<double> init_param) 
{
    return(1.);
}

double initial_null(vector<double> x,vector<double> init_param) 
{
    double a = 0.;
    return(a);
}



double initial_gauss(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( init_param*exp(-pow((x*dev_std),2)) );
}

double initial_gauss_2D(vector<double> x,vector<double> init_param)
{
    double amplitude = init_param[0];
    double dev_std = 1;//init_param[1];
    return( amplitude*exp(-pow((x[0]*dev_std),2)) * -pow((x[1]*dev_std),2)) );
}

double initial_gauss_PHI(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( -2.*x*dev_std*dev_std*init_param*exp(-pow((x*dev_std),2)) );
}

double initial_gauss_PHI_compactified(vector<double> x,vector<double> init_param)
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

double initial_gauss_phi_compactified(vector<double> x,vector<double> init_param)
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


// ------ hyperboloidal compactification and Chi function rescaling  ------ //

double initial_gauss_PI_compactified_Chi1(double r,vector<double> init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = -dev_std*dev_std*(-r+2.*R)*(-r+2.*R);
    double C = dev_std*dev_std*r*r;
    double D = (-r+2.*R)*(-r+2.*R);
    if(r!=0 && r!=s)
    {
        return( (-A*exp(-C)+A*exp(B)+2.*A*C*exp(-C)-2.*A*D*ds_sq*exp(B))*Chi/2/R);
    }
    else if (r==0)
    {
        return(0);
    }
    else
    {
        return(1./2.*(-A*exp(-C)+2*A*C*exp(-C)));
    }
}
double initial_gauss_PHI_compactified_Chi(vector<double> x,vector<double> init_param)
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

double initial_gauss_PHI_compactified_Chi1(double r,vector<double> init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = exp(-pow(dev_std,2)*pow(r-2*R,2));
    double C = r-2*R;
    double D = -r+R;
    double G = dev_std*dev_std*r*r;
    if(r!=s && r!=0)
    {
        return( (A*(exp(-G)*(D-2*G*R)+B*(r-R*(1+2*G+8*D*dev_std*dev_std*R)))*Chi+(-A*B*C+A*exp(-G)*r)*R*ChiPrime)/(2*R*R) );
    }    
    else if(r==0)
    {
        return(0);
    }
    else
    {
        return(1./2.*exp(-pow(dev_std*s,2))*(A-2.*A*pow(dev_std*s,2)));
    }
}

double initial_gauss_phi_compactified_Chi(double r,vector<double> init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);

    if(r!=s)
    {
        return( init_param*Chi*exp(-pow(dev_std*R,2)));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_phi_compactified_Chi1(double r,vector<double> init_param)
{
    double s = 5.;
    double dev_std = 1./5.;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);

    if(r!=s && r!=0)
    {
        return((init_param*exp(-pow(dev_std*r,2))*r + init_param*exp(-pow(dev_std*(-r+2.*R),2))*(-r+2.*R))/2./R * Chi);
    }
    
    else if(r==0)
    {
        return(1);
    }
    else
    {
        return(1./2.* exp(-pow(dev_std*s,2))*s);
    }
}


// ------  ------  ------  ------  ------  ------  ------  ------  ------  ------


double initial_parab_function(vector<double> x,vector<double> init_param)
{
    double a = -2*pow(x,2);
    return(a);
}

double line(vector<double> x,double t)
{
    return(4*t);
}



// initialize the grid function
grid initialize_grid(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim)
{
    
    // setting the ghost point equal to the order of the artificial dissipator
    if (ord>gl)
    {
        gl = ord;
    }
    if(ord>gr)
    {
        gr = ord;
    }
    
    // initializing the size of all the uni-dimensional grid
    
    vector<int> i_v(dim);
    for(int d=0;d<dim;d++)
    {
        if(dmax[d]!=dmin[d])
        {
            i_v[d] =  int( ((dmax[d]+dx[d]*double(gr))-(dmin[d]-dx[d]*double(gl)) +dx[d]/2 ) / dx[d]) + 1;
        }
        else
        {
            i_v[d] = 1;
        }
            cout<<i_v[d]<<endl;
    }

    grid Grid;

    double val1;
    double val2;
    double val3;
    for(int i=0;i<i_v[0];i++)
    {
        vector<vector<vector<double>>> support1 ;
        Grid.push_back(support1);
        for(int j=0;j<i_v[1];j++)
        {
            vector<vector<double>> support2 ;
            Grid[i].push_back(support2);
            for(int k=0;k<i_v[2];k++)
            {
                val1 = dmin[0]+dx[0]*(i-gl);
                val2 = dmin[1]+dx[1]*(j-gl);
                val3 = dmin[2]+dx[2]*(k-gl);
                Grid[i][j].push_back({val1, val2, val3});
            }
        }
    }
    return(Grid);
} 


// initialize the scalar fields
fields_vector initialize_fields(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim,grid Grid,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> &param_ic)
{
    int N = funcs.size();
    fields_vector new_fields;
    // setting the ghost point equal to the order of the artificial dissipator
    if (ord>gl)
    {
        gl = ord;
    }
    if(ord>gr)
    {
        gr = ord;
    }
    
    for(int n=0;n<N;n++)
    {
        vector<vector<vector<double>>> support1 ;
        new_fields.push_back(support1);
        for(int i=0;i<Grid.size();i++)
        {
            vector<vector<double>> support2 ;
            new_fields[n].push_back(support2);
            for(int j=0;j<Grid[i].size();j++)
            {
                vector<double> support3;
                new_fields[n][i].push_back(support3);
                for(int k=0;k<Grid[i][j].size();k++)
                {
                    cout<<"flag";
                    new_fields[n][i][j].push_back(funcs[n](Grid[i][j][k],param_ic));
                }
            }
        }
    }
    return(new_fields);
} 
/*
// initialize the fields at time zero
vector<vector<double>> initialize_fields(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> &param_ic,int gl, int gr,int ord,int dim)
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
    vector<int> S(dim);
    for(int d=0;d<dim;d++)
    {
        S[d]=  int( ((dmax[d]+dx[d]*double(gr))-(dmin[d]-dx*double(gl)) +dx[d]/2 ) / dx[d]) + 1;
    }
    //cout<<S<<endl;
    std::vector<std::vector<double>> new_vect(N,dim, vector<double>(S));
    for(int n=0;n<N;n++)
    {
        for(int d=0;d<dim;d++)
        {
            new_vect[n][d].push_back(vector<double>(S[d]));
        }
    }
    
    
    double val;
    
    for(int d=0;d<dim;d++)
    {
        for(double j=gl;j < S-gr ;j=j+1)
        {
            val = d_min[d]+dx[d]*(j-gl);
            //cout<< j <<"\n";
            new_vect[0][d][j] = val
        }
    }
    
    for(int d=0;d<dim;d++)
    {
        for(int n=0;n<N;n++)
        {
            for(double j=gl;j < S-gr ;j=j+1)
            {
                val = d_min[d]+dx[d]*(j-gl);
                //cout<< j <<"\n";
                new_vect[n][d][j] = funcs[n](val,param_ic);
            }
        }
    }
    return(new_vect);
}
*/
/*
void init_func(std::vector<double> &func_vect,vector<double> &dmin,vector<double> &dmax,vector<double> &dx,double(*func)(vector<double>,vector<double>),double t)
{
    for(int d=0;d<dim;d++)
    {
        for(double j=dmin[d];j<dmax[d]+dx/2.;j+=dx[d])
        {
            func_vect.push_back(func(j,t));
        }
    }
}
*/
// --------- BOUNDARY CONDITIONS --------- //

// adv eq 
void adv_boundaries_TEM(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //same of the bulk
    fields_vect_new[j].insert(fields_vect_new[j].begin(),((-4*fields_vect_old[0][gl]+7.*fields_vect_old[0][gl+1]-4.*fields_vect_old[0][gl+2]+fields_vect_old[0][gl+3])/2/dx));
    //boundary conditions
    fields_vect_new[j].push_back((-(4*fields_vect_old[0][last_ind]-7.*fields_vect_old[0][last_ind-1]+4.*fields_vect_old[0][last_ind-2]-fields_vect_old[0][last_ind-3]) ));
    //fields_vect_new[j].push_back(((5.*fields_vect_old[j][last_ind]-11.*fields_vect_old[j][last_ind-1]+10.*fields_vect_old[j][last_ind-2]-5*    fields_vect_old[j][last_ind-3]+fields_vect_old[j][last_ind-4])/2/dx ));
}

void adv_boundaries_right(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //same of the bulk
    fields_vect_new[j][gl] = Dx[0](fields_vect_old[0],gl,dx) + artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt);
    //boundary conditions
    fields_vect_new[j][last_ind] = (-param[0]*Dx[0](fields_vect_old[0],last_ind,dx));
}

void adv_boundaries_left(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx)+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    //fields_vect_new[j].insert(fields_vect_new[j].begin(),(-3./2.*fields_vect_old[0][gl]+2.*fields_vect_old[0][gl+1]-1*2.*fields_vect_old[0][gl+2])+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    
    

    // same of the bulk
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx));
    //fields_vect_new[j].push_back((-(3./2.*fields_vect_old[0][last_ind]-2.*fields_vect_old[0][last_ind-1]+1/2.*fields_vect_old[0][last_ind-2] )));
}



void adv_boundaries_periodic(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //left boundary
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl+1]-fields_vect_old[0][last_ind])/2/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind-1])/2/dx ));
    
}

void adv_boundaries_periodic_backward_der(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    //left boundary
    
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind])/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));   
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][last_ind]-fields_vect_old[0][last_ind-1])/dx ));

}


// wave equation

void refl_abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( (Dx[0](fields_vect_old[0],gl,dx)-Dx[0](fields_vect_old[1],gl,dx))/2 ));
    fields_vect_new[j].push_back((-(Dx[0](fields_vect_old[0],last_ind,dx)-Dx[0](fields_vect_old[1],last_ind,dx))/2));
}



void radiative_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    int ind_space = gl;
    // left boundary conditions
    fields_vect_new[j][gl] = ( Dx[0](fields_vect_old[0],gl,dx));

    //fields_vect_new[j].insert(fields_vect_new[j].begin(), Dx[0](fields_vect_old[0],gl,dx));
    //fields_vect_new[j].insert(fields_vect_new[j].begin(),( (-3/2*fields_vect_old[0][gl]+2*fields_vect_old[0][gl+1]-1/2*fields_vect_old[0][gl+2])/dx +artificial_diss(epsilon,ord,fields_vect_old,0,gl,dx,dt)));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx));

    //fields_vect_new[j].push_back(-Dx[0](fields_vect_old[0],last_ind,dx));
    //fields_vect_new[j].push_back(- ((3/2*fields_vect_old[0][last_ind]-2*fields_vect_old[0][last_ind-1]+1/2*fields_vect_old[0][last_ind-2])/dx +artificial_diss(epsilon,ord,fields_vect_old,0,last_ind,dx,dt)));
}

void radiative_outer_boundaries_PI_we(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
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


void abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( Dx[0](fields_vect_old[0],gl,dx)));
    fields_vect_new[j].push_back((-Dx[0](fields_vect_old[0],last_ind,dx)));
}

void abs_boundaries_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j].insert(fields_vect_new[j].begin(), (Dx[0](fields_vect_old[1],gl,dx)));
    // right boundary conditions
    fields_vect_new[j].push_back(( -Dx[0](fields_vect_old[1],last_ind,dx)));
}

   
void no_boundary_conditions_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[0],last_ind,dx));
}

void no_boundary_conditions_phi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[0][gl]);
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[0][last_ind]);
}

// compactified and rescaled wave equation//

// R rescaling
void no_boundary_conditions_PI_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
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
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);    
}

// R rescaling
void no_boundary_conditions_PHI_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
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
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);    
}

// R rescaling
void no_boundary_conditions_phi_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    
    double left_value = -fields_vect_old[0][ind_space_left];
    double right_value = -fields_vect_old[0][last_ind];
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = right_value;    
}


// CHI RESCALING //


void no_boundary_conditions_PI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    /*
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
    */
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    /*
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
    */
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_PHI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    /*
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
    */
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    /*
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
     */
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_phi_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
        
    double r = dmin;   
    
    //double left_value = fields_vect_old[0][ind_space_left];
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    r = dmax;

    //double right_value = fields_vect_old[0][last_ind];
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}
// ------------- // MODEL 1 // ------------- //


void radiative_outer_boundaries_PI_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    
    
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3.* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +(pow(fields_vect_old[1][ind_space],2)-pow(fields_vect_old[0][ind_space],2) );
    
    fields_vect_new[j][ind_space] = ( derivative);

    
    
    // right boundary conditions
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx));

    
}



// ------------- // MODEL 3 // ------------- //


void radiative_outer_boundaries_PI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))      
                        + (fields_vect_old[2][ind_space] + param[0] * fields_vect_old[5][ind_space])/pow(param[0],2)*
                        (pow(fields_vect_old[1][ind_space],2)+pow(fields_vect_old[4][ind_space],2)-pow(fields_vect_old[0][ind_space],2)-pow(fields_vect_old[3][ind_space],2));
    
    fields_vect_new[j][ind_space] = ( derivative);

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)-fields_vect_old[0][last_ind]/x );
}

void radiative_outer_boundaries_PI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = 3* ( pow((x+dx),2)*fields_vect_old[4][ind_space+1] - pow((x-dx),2)*fields_vect_old[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) 
                        + (fields_vect_old[5][ind_space] - param[0] * fields_vect_old[2][ind_space])/pow(param[0],2)*
                        (pow(fields_vect_old[1][ind_space],2)+pow(fields_vect_old[4][ind_space],2)-pow(fields_vect_old[0][ind_space],2)- 
                        pow(fields_vect_old[3][ind_space],2));
    
                        
                        
    fields_vect_new[j][ind_space] = ( derivative);

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[3],last_ind,dx)-fields_vect_old[3][last_ind]/x );
}


void no_boundary_conditions_PHI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[0],last_ind,dx));
}

void no_boundary_conditions_PHI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[3],gl,dx));
    
    // right boundary conditions
    fields_vect_new[j][last_ind] = ( Dx[0](fields_vect_old[3],last_ind,dx));
}


void no_boundary_conditions_phi1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[0][gl]);
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[0][last_ind]);
}

void no_boundary_conditions_phi2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    fields_vect_new[j][gl] = (fields_vect_old[3][gl]);
    // right boundary conditions
    fields_vect_new[j][last_ind] = (fields_vect_old[3][last_ind]);
}

//-------------- DIFFERENTIAL OPERATORS -------------//
    
// simple function to get the derivative of a function, f=g'
// It uses boundary condition!
// input:primitive function (g),discretization step


    
// fourth order central finite difference expression for a first spatial derivative
double first_der_fourth_order_centered(std::vector<double> vector,int i,vector<double> dx)
{
    return((vector[i-2]/12.-2./3.*vector[i-1]+2./3.*vector[i+1]-vector[i+2]/12.)/dx);
}
    
// fourth order central finite difference expression for a first spatial derivative
double first_der_second_order_centered(std::vector<double> vector,int i,vector<double> dx)
{
    return((-vector[i-1]+vector[i+1])/2./dx);
}

double first_der_first_order_backward(std::vector<double> vector,int i,vector<double> dx)
{
    return((-vector[i-1]+vector[i])/dx);
}

double first_der_second_order_backward(std::vector<double> vector,int i,vector<double> dx)
{
    return((1/2.*vector[i-2]-2.*vector[i-1]+3./2.*vector[i])/dx);
}

double first_der_second_order_forward(std::vector<double> vector,int i,vector<double> dx)
{
    return((1/2.*vector[i-2]+2.*vector[i-1]-3./2.*vector[i])/dx);
}
    
//-------------- GHOST POINTS --------------//

void ghost_point_extrapolation_1_ord(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation1(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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


void ghost_point_extrapolation2_TEM(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_2_ord(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
{
    int dim = 3
    // attachment of the ghost points to the boundary of the function
    for(int d;d<dim;d++)
    {
        for (int i=gl-1; i>=0;i--)
        {
            
            field_vect[j][i][j][k] = field_vect[j][d][i+1]-(-25./12.*field_vect[j][d][i+1]+4.*field_vect[j][d][i+2]-3.*field_vect[j][d][i+3]+4./3.*field_vect[j][d][i+4]-1./4.*field_vect[j][d][i+5])+1./2.*(2.*field_vect[j][d][i+1]-5.*field_vect[j][d][i+2]+4.*field_vect[j][d][i+3]-field_vect[j][d][i+4]);
        }
        
        for (int i=field_vect[j][d].size()-gr; i<field_vect[j][d].size();i++)
        {
            field_vect[j][d][i] = field_vect[j][d][i-1]+(+25./12.*field_vect[j][d][i-1]-4.*field_vect[j][d][i-2]+3.*field_vect[j][d][i-3]-4./3.*field_vect[j][d][i-4]+1./4.*field_vect[j][d][i-5]+1./2.*(2.*field_vect[j][d][i-1]-5.*field_vect[j][d][i-2]+4.*field_vect[j][d][i-3]-field_vect[j][d][i-4]));
        }
    }
}

void ghost_point_extrapolation_6_ord(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation3(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_1_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_4_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_1_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_2_ord_spline_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
{
    int S = field_vect[j].size();
    // computation of the ghost points via interpolation
    // we find the x vector and the relative y vector values of the subsections starting from the boundaries
    int number_of_point_int = 5;        

    std::vector<double> right_side_x;
    std::vector<double> right_side_y;
    for (int i = 0;i<number_of_point_int+1;i++)
    {
        right_side_x.insert(right_side_x.begin(),(dmax-i*dx));
        //cout<<"x:"<<(dmax-i*dx);
        right_side_y.insert(right_side_y.begin(),(field_vect[j][(S-1)-i-gr]));
        //cout<<" y:"<<(field_vect[j][field_vect[j].size()-1-i-gr]);
    }
    
    std::vector<double> left_side_x;
    std::vector<double> left_side_y;
     for (int i = 0;i<number_of_point_int+1;i++)
    {
        left_side_x.push_back(dmin+i*dx);
        //cout<<"x:"<<(dmax-i*dx);
        left_side_y.push_back(field_vect[j][gl+i]);
        //cout<<" y:"<<(field_vect[j][field_vect[j].size()-1-i-gr]);
    }
    
    // interpolation 
    tk::spline right(right_side_x,right_side_y);
    tk::spline left(left_side_x,left_side_y);
    
    //field_vect[j][gl] = left(dmin);
    
    
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
    
    int counter = 1;
    for (int i=S-gr; i<S;i++)
    {
        
        field_vect[j][i] = right(dmax+dx*counter);
        counter += 1;
    }
    /*
    for (int i=0; i<gl;i++)
    {
        
        field_vect[j][i] = left(dmin-dx*gl+dx*i);
    }*/
    
}

void ghost_point_extrapolation_4_ord_spherical_symmetry_rescaled(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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

void ghost_point_extrapolation_6_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax)
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
        field_vect[j][i]=(field_vect[j][i]+(49./20.*field_vect[j][i]-6.*field_vect[j][i-1]+15./2.*field_vect[j][i-2]-20./3.*field_vect[j][i-3]+15./4.*field_vect[j][i-4]-6./5.*field_vect[j][i-5]+1./6.*field_vect[j][i-6]));
    }
}


// ------------- ARTIFICIAL DISSIPATION -------------- //
double artificial_dissipation_2(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,vector<double> dx,double dt)      
{
    return(-epsilon*(dt/dx)*pow(-1,ord)*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_Husa(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,vector<double> dx,double dt)      
{
    return(-epsilon/pow(dx,1)*pow(-1,ord)/16*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
}

// to control!
double artificial_dissipation_4(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,vector<double> dx,double dt)      
{
    return(-epsilon*pow(dx,3)*pow(-1,ord)/4*(copy_fields_vect[j][i+2]-4*copy_fields_vect[j][i+1]+6*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_ani(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,vector<double> dx,double dt)      
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

double conv_test(std::vector<double> &num1,std::vector<double> &num2,double(*theo_sol)(double,double),void(*init_func)(std::vector<double> &,double ,double ,double ,double(*)(double,double),double),double h1,double h2, double(*norm) (std::vector<double> & , double ),double (*norm_of_diff)(std::vector<double> &, std::vector<double> &,double , double(*) (std::vector<double> & , double )),vector<double> &dmin, vector<double> &dmax,vector<double> dx,double t)
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



void print_f(fields_vector fields_vect, vector<double> &dmin, vector<double> &dmax, vector<double> dx, string name_file,string name_folder, int gl, int gr,MPI_Status status, int totalnodes, int mynode,MPI_Request request)
{
    int dim = 3;
    vector<int> S(dim);
    vector<double> index_dmin_local(dim),index_dmax_local(dim);
    for(int d=0;d<dim;d++)
    {
        S[d] =  int( ((dmax+dx*double(gr))-(dmin-dx*double(gl)) +dx/2 ) / dx) + 1;
        // defining the first and last spatial index of the domain, the actual process will work inside this range
            
        index_dmin_local[d] = (gl+1) + mynode * int((S[d]-2-gl-gr)/totalnodes);
        if (mynode==totalnodes-1)
        {
            index_dmax_local[d] =  S[d]-1-gr;
        }
        else
        {
            index_dmax_local[d] = (gl+1) + (mynode+1) * int((S[d]-2-gl-gr)/totalnodes);
        }
    }
    
    
    ofstream myfile_print_f;
    myfile_print_f.open (name_file,ios::app);
    
    // headers of the columns
    for(int d=0;d<dim;d++)
    {
        myfile_print_f<<"x_"<<to_string(d)<<",";
    }
    for (int i=0; i<fields_vect.size()-1; i++)
    {
        myfile_print_f << "field"<<to_string(i)<<",";
    }    
    myfile_print_f << "field"<<to_string(fields_vect.size()-1)<<"\n";
    
    
    // Each processor print only the spatial points that it is evolving
    
    // the processor zero print also the left ghost points
    if(mynode==0)
    {
        for(int j=0;j<gl+1;j++)
        {
            for(int d=0;d<dim;d++)
            {
                myfile_print_f << dmin[d]+dx[d]*double(j-gl)<<",";
            // for every fields add the relative value
            for (int i=0; i<fields_vect.size()-1; i++)
            {
                myfile_print_f << fields_vect[i][j]<<",";
            }
            myfile_print_f << fields_vect[fields_vect.size()-1][j];
            myfile_print_f<<"\n";
            // send a message from proc 0 to the other, saying that the header has been written
            
        }
    }
    
    
    
    for (int j=index_dmin_local; j<index_dmax_local; j= j+1)
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
    
    // the last processor print also the right ghost points
    if(mynode==totalnodes-1)
    {
        for(int j=index_dmax_local;j<fields_vect[0].size();j++)
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
    }
    
    myfile_print_f.close();
}


void read_parameters(string name_parameters_file, vector<double> &dmin, vector<double> &dmax, vector<double> &h1, double &integration_interval,int &step_to_save,int &gl,int &gr,int &ord, vector<double> &epsilon,vector<double> &parameters_ic_vector, vector<double> &parameters)
{
    // we read the input parameter from an external file
    ifstream input_file;
    input_file.open(name_parameters_file);
    
    // dmin vector
    std::string dmin_v;
    input_file.ignore(256,' ');
    getline (input_file,dmin_v);
    
    std::stringstream iss3 (dmin_v);
    double number;
    while ( iss3 >> number )
    {
        dmin.push_back( number );
    }
    
    
    // dmin vector
    std::string dmax_v;
    input_file.ignore(256,' ');
    getline (input_file,dmax_v);
    
    std::stringstream iss4 (dmax_v);
    while ( iss4 >> number )
    {
        dmax.push_back( number );
    }
    
    // dmin vector
    std::string h1_v;
    input_file.ignore(256,' ');
    getline (input_file,h1_v);
    
    std::stringstream iss5 (h1_v);
    while ( iss5 >> number )
    {
        h1.push_back( number );
    }
    
    
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
    //double number;
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


// Function usefull in a parallelized setting

void communication(std::vector< std::vector<double> > &fields_vect,int j,int index_dmax_local,int index_dmin_local,int nitems,int mynode,int totalnodes,MPI_Status status)
{
    if(nitems==1)
    {
    // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
    // each processor needs to receive (and send) the borders of its subdomain
        if(mynode==0)
        {
            MPI_Sendrecv(&(fields_vect[j][index_dmax_local-1]),nitems, MPI_DOUBLE, mynode+1,
            0, &(fields_vect[j][index_dmax_local]), nitems, MPI_DOUBLE, mynode+1,  0, MPI_COMM_WORLD, &status);
        }
            
        if(mynode==totalnodes-1)
        {
            MPI_Sendrecv(&(fields_vect[j][index_dmin_local]),nitems, MPI_DOUBLE, mynode-1,
            0, &(fields_vect[j][index_dmin_local-1]), nitems, MPI_DOUBLE, mynode-1,  0, MPI_COMM_WORLD, &status);
        }
            
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            MPI_Send(&fields_vect[j][index_dmax_local-1],nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD);
            MPI_Send(&fields_vect[j][index_dmin_local],nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD);
            MPI_Recv(&(fields_vect[j][index_dmax_local]),nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD, &status);
            MPI_Recv(&(fields_vect[j][index_dmin_local-1]),nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD, &status);
        }
    }
    else
    {
        
        if(mynode==0)
        {
            double message_out[nitems] = {fields_vect[j][index_dmax_local-1],fields_vect[j][index_dmax_local-2]};
            double message_in[nitems] ;
            MPI_Sendrecv(&message_out,nitems, MPI_DOUBLE, mynode+1,
            0, &message_in, nitems, MPI_DOUBLE, mynode+1,  0, MPI_COMM_WORLD, &status);
            for(int i=0;i<nitems;i++)
            {
            fields_vect[j][index_dmax_local+i] = message_in[i];
            }
        }
        else if(mynode==totalnodes-1)
        {
            double message_out[nitems] = {fields_vect[j][index_dmin_local],fields_vect[j][index_dmin_local+1]};
            double message_in[nitems];
            MPI_Sendrecv(&message_out,nitems, MPI_DOUBLE, mynode-1,
            0, &message_in, nitems, MPI_DOUBLE, mynode-1,  0, MPI_COMM_WORLD, &status);
            for(int i=0;i<nitems;i++)
            {
            fields_vect[j][index_dmin_local-1-i] = message_in[i];
            }
        }
        else
        {
            double message_out_left[nitems] = {fields_vect[j][index_dmin_local],fields_vect[j][index_dmin_local+1]};
            MPI_Send(&message_out_left,nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD);
            
            double message_out_right[nitems] = {fields_vect[j][index_dmax_local-1],fields_vect[j][index_dmax_local-2]};
            MPI_Send(&message_out_right,nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD);
            
            double message_in_left[nitems];
            MPI_Recv(&message_in_left,nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD, &status);
            for(int i=0;i<nitems;i++)
            {
            fields_vect[j][index_dmin_local-1-i] = message_in_left[i];
            }
            
            double message_in_right[nitems];
            MPI_Recv(&message_in_right,nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD, &status);
            for(int i=0;i<nitems;i++)
            {
            fields_vect[j][index_dmax_local+i] = message_in_right[i];
            }
        }
    }
}
