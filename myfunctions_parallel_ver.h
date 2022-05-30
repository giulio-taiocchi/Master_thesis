#include <math.h>
#include <fstream>
#include <sstream>
#include <mpi.h>
#include "spline.h"
using namespace std;

typedef void (*print_function)(std::vector< std::vector<double> > fields_vect, double dmin,double dmax, double dx, string name_file,string name_folder, int gl, int gr,MPI_Status status, int totalnodes, int mynode,MPI_Request request);


typedef double (*artificial_dissipation_function)(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt) ;

typedef std::vector<double (*)(std::vector<double>,int,double)>  derivative_vector;

typedef void(*ghost_point_extrapolation_function)(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param);

typedef double (*evolution_function)(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl);

typedef void(*boundary_conditions_function)(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo);

typedef void(*communication_function)(std::vector< std::vector<double> > &fields_vect,int j,int index_dmax_local,int index_dmin_local,int nitems,int mynode,int totalnodes,MPI_Status status);

typedef vector<vector<double>>(*one_step_function)(std::vector< std::vector<double> > fields_vect,double dmin,double dmax,double dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord,MPI_Status status, int totalnodes, int mynode, communication_function communication)   ;

typedef void(*method_of_line_function)(std::vector< std::vector<double> > fields_vect,one_step_function one_step, double dx, std::vector<double> param, double dt, double interval,    double dmin,    double dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path ,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication);

typedef vector<vector<double>>(*initialization_function)(double d_min,double d_max,double dx,std::vector<double(*)(double,double,vector<double>)> funcs,double param_ic,vector<double> param,int gl, int gr,int ord);

void multiple_parameters_run(std::vector<double>& parameters_ic_vector, std::vector<double(*)(double, double,vector<double>)>& initial_conditions, initialization_function initialize_fields, double dmin, double dmax, double h1, double h2, double h3, double dt1, double dt2, double dt3, double integration_interval,int step_to_save, derivative_vector & Dx ,std::vector< evolution_function > & R_vector, std::vector< boundary_conditions_function > & b_func,  std::vector<double> & parameters, one_step_function one_step, int gl,int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2, vector<double> epsilon1,print_function print_f, string file_path, method_of_line_function MOL_RK4,int ord,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication)
{
    ofstream myfile2;
    //myfile2.open (file_path+"name_of_file");
    //myfile2<<"names\n";
    //myfile2.close();
    
    for(int e=0;e<epsilon1.size();e++)
    {
        for (int l=0;l<parameters_ic_vector.size();l++)
        {
            if(mynode==0)
            {
                cout<<"------------------------------------------------------------------------------------------"<<endl;
                cout<<"processor :"<<mynode<<"\nepsilon = "<<epsilon1[e]<<" amplitude = "<<parameters_ic_vector[l]<<endl;
                cout<<"------------------------------------------------------------------------------------------"<<endl;
            }
            // initialize fields
            
            
            
            std::vector< std::vector<double> > fields_vect1 = initialize_fields(dmin,dmax,h1,initial_conditions,parameters_ic_vector[l],parameters,gl,gr,ord);     
            std::vector< std::vector<double> > fields_vect2 = initialize_fields(dmin,dmax,h2,initial_conditions,parameters_ic_vector[l],parameters,gl,gr,ord);
            std::vector< std::vector<double> > fields_vect3 = initialize_fields(dmin,dmax,h3,initial_conditions,parameters_ic_vector[l],parameters,gl,gr,ord);
            // setting the output variables
            //cout<<"\n lun vec \n"<<endl;
            string name_file = "processor_"+to_string(mynode)+"_ampl_"+to_string(parameters_ic_vector[l])+"_eps"+to_string(epsilon1[e])+"_dx_"+to_string(h1)+"steps"+to_string(step_to_save)+"last_time"+to_string(integration_interval)+".csv";
            
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
    

void MOL_RK4(std::vector< std::vector<double> > fields_vect,one_step_function one_step, double dx, std::vector<double> param, double dt, double interval,    double dmin,    double dmax,std::vector< evolution_function > R_vect,std::vector< boundary_conditions_function > bc, double step_to_save,print_function print_f,int gl, int gr,ghost_point_extrapolation_function ghost_point_extrapolation,artificial_dissipation_function artificial_diss_2,double epsilon,int ord,derivative_vector Dx,string file_path,MPI_Status status, int totalnodes, int mynode,MPI_Request request, communication_function communication)
{   
    cout<<"\nProcessor :"<<mynode<<endl<<"--- Method of lines called ---\ndx = "<<dx<<"\ndt = "<<dt<<"\nDomain = ["<<dmin<<","<<dmax<<"]\nlast time objective :"<<interval<<"\n";
    
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
        ghost_point_extrapolation(fields_vect, 0,dx,dt,j,gl,gr,dmin,dmax,param);
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


vector<vector<double>> onestep_RK4_1(std::vector< std::vector<double> > fields_vect,double dmin,double dmax,double dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,std::vector< boundary_conditions_function > bc,double t, int gl, int gr, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon,derivative_vector Dx,int ord,MPI_Status status, int totalnodes, int mynode, communication_function communication)  
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
    
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,param);
    }
    
    // k1 building
    for (int j=0; j <N; j++)
    {
        // evualuating the "bulk" of k1
        // we have to consider the physical domain of the fields (so exclude the GP) and then exclude the boundaries value (+1 and -1)
        // we are also dividing the spatial domain in subdomains, in which each processor run in parallel
        
        // calculating the k1 function (just in the interval associated to the actual processor)
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k1[j][i] = (evo[j](j,i,fields_vect,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl))
                            +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        communication(k1,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
    
        
        // evaluating the boundary of k1
        bc[j](k1,fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        k1[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt);
        k1[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt);
        // GP extrapolation for k1 vector
    }
    for (int j=0; j <N; j++)
    {
    // computing the argument for the next coefficient k2
        if(mynode==0)
        {
            ghost_point_extrapolation(k1, t,dx,dt,j,gl,gr,dmin,dmax,param);
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }
        }
        if(mynode==totalnodes-1)
        {
            ghost_point_extrapolation(k1, t,dx,dt,j,gl,gr,dmin,dmax,param);
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }
        }
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k1[j][i] = (k1[j][i])*dt/2. + fields_vect[j][i];//+ artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }
        }
    }
    
    
    
    // k2 building
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k2[j][i] = evo[j](j,i,support_k1,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt/2,gl)
            +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(k2,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
        //boundary conditions update
        bc[j](k2,support_k1, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        // artificial dissipation
        k2[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt/2);
        k2[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt/2);
        
    }
    for (int j=0; j <N; j++)
    {
        
        if(mynode==0)
        {
            ghost_point_extrapolation(k2, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,param);
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];//+ artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }  
        }
        if(mynode==totalnodes-1)
        {
            ghost_point_extrapolation(k2, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,param);
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }
        }
        
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k2[j][i] = (k2[j][i])*dt/2. + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2);
            }
        }
    }
    
    // k3 building
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k3[j][i] = evo[j](j,i,support_k2,dx,dmin,param,t+dt/2.,Dx,artificial_diss,epsilon,ord,dt,gl)
            +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(k3,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        bc[j](k3,support_k2, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        k3[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt/2);
        k3[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt/2);
    }
    for (int j=0; j <N; j++)
    {
        
        if(mynode==0)
        {
            ghost_point_extrapolation(k3, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,param);
            for (int i=0;i<index_dmax_local+1;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
        }
        
        if(mynode==totalnodes-1)
        {
            ghost_point_extrapolation(k3, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,param);
            for (int i=index_dmin_local-1;i<S;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];// + artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
        }
        
        if(mynode!=0 && mynode!=totalnodes)
        {
            for (int i=index_dmin_local-1;i<index_dmax_local+1;i++)
            {
                support_k3[j][i] = k3[j][i]*dt + fields_vect[j][i];//+ artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
            }
        }
        
    }
   
   
    // k4 building    
    for (int j=0; j <N; j++)
    {
        for (int i=index_dmin_local;i<index_dmax_local;i++)
        {
            k4[j][i] = evo[j](j,i,support_k3,dx,dmin,param,t+dt,Dx,artificial_diss,epsilon,ord,dt/2,gl)
            +artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt);
        }
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        communication(k4,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        
        bc[j](k4,support_k3, t+dt,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        k4[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt);
        k4[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt);
    }
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(k4, t+dt,dx,dt,j,gl,gr,dmin,dmax,param);
    }
    // we create a new vector that contains all the new fields. It is a support vector that will be swapped with the old one
    
    std::vector<std::vector<double>> new_fields_vect(N, vector<double>(S));

    for (int j=0; j <N; j++)
    {
        if(mynode==0)
        {
             for (int i=gl;i<index_dmax_local;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2.*k2[j][i]+2.*k3[j][i]+k4[j][i])/6.;
                //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/6);
            }               
        }
        
        if(mynode==totalnodes-1)
        {
            for (int i=index_dmin_local;i<S;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2.*k2[j][i]+2.*k3[j][i]+k4[j][i])/6.;
                                        //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/6);
            }            
        }
        
        if(mynode!=0 && mynode!=totalnodes-1)
        {
            for (int i=index_dmin_local;i<index_dmax_local;i++)
            {
                new_fields_vect[j][i] = fields_vect[j][i] + dt*(k1[j][i]+2.*k2[j][i]+2.*k3[j][i]+k4[j][i])/6.;
                                        //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/6);
            }
        }
        communication(new_fields_vect,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
    }   
    

    // populate the ghost zone of the fields, we need them to calculate k1 in the next time step
    for (int j=0; j <N; j++)
    {
        if(mynode==0 || mynode==totalnodes-1)
        {
            ghost_point_extrapolation(new_fields_vect, t,dx,dt,j,gl,gr,dmin,dmax,param);
        }
    }
    //cout<<"old "<<fields_vect.size()<<"new "<<new_fields_vect.size()<<endl;
    
    return(new_fields_vect);
}



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
    return (3*(pow((x+dx),2)*fields_vect[1][ind_space+1]-pow((x-dx),2)*fields_vect[1][ind_space-1])/(pow((x+dx),3)-pow((x-dx),3)) );
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


//------------------- WAVE EQUATION: CHI RESCALING HYPERBOLOIDAL FOLIATION -------------------//

// compactification: T = t + H(R), H'=1-1/R', R=r/(1-(r/s)^2)
// rescaling: by chi(R) = (1+R^2)^(1/2)
// note: at the oriign we apply Evan's method for the terms wich present 1/r
// note: for r=s, we use a different RHS sicne some of the terms would diverge. Computing a limit "by hand", we can put = 0 the problematic terms

double wave_eq_compactified_PI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2. + s*s);
    double B = pow(r,4) - 4. * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3.*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3.*s*s;
    double H = -pow(r,4) + 4.*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1. + s*s);
    double O = pow(r,4)-pow(s,4);
    
    return (2.*M*r*s*s*fields_vect[1][ind_space]/A/B
            +3.*O*O*pow(s,4)*fields_vect[2][ind_space]/(pow(A,2)*B)
            +D*Dx[0](fields_vect[0],ind_space,dx)/H
            -N*r*r*s*s*Dx[0](fields_vect[1],ind_space,dx)/A/B
            -pow(s,8)/A/B* 3.* (pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3))
            );
}

double wave_eq_compactified_PHI_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*(ind_space-gl);
    double s = param[0];
    
    double A = pow(r,4) + pow(s,4) + r*r*s*s *(-2. + s*s);
    double B = pow(r,4) - 4. * pow(r,2)* s*s - pow(s,4);
    double C = r- s;
    double D = pow(r,4) - 3.*r*r*s*s;
    double F = pow(r,2) + s*s;
    double G = pow(r,2)-3.*s*s;
    double H = -pow(r,4) + 4.*pow(r,2)*s*s + pow(s,4);
    double L = pow(r,3)*s - r* pow(s,3);
    double M = pow(r,4) + r*r * s*s - pow(s,4);
    double N = pow(r,4) - pow(s,4) + pow(s,6) + r*r* s*s* (-1. + s*s);
    double O = pow(r,4)-pow(s,4);
    
    return ( 2.*C*F*G*r*(r+s)*fields_vect[1][ind_space] /(A*B)
            +(3.*F*G*L*L*fields_vect[2][ind_space] )/(A*A*B)
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

// ----------- // MODEL 1Bis // ----------- //

double model1_PI_bis(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3.* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +param[0]*pow(x*x+1.,0.5)*(pow(fields_vect[1][ind_space],2)-pow(fields_vect[0][ind_space],2)) );
}

//------------------- MODEL 1: CHARACTERISTIC VARIABLE IN FLAT SPACE -------------------//
//fields_vect [0]-> PsiPlus, [1]->PsiMinus, [2]->Psi

double model1_charvar_Psi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    return((fields_vect[0][ind_space]+fields_vect[1][ind_space])/2.);
}

double model1_charvar_PsiPlus(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    return(-fields_vect[0][ind_space]*fields_vect[1][ind_space]
            +Dx[0](fields_vect[0],ind_space,dx)/2.
            +Dx[0](fields_vect[1],ind_space,dx)/2.
            +0.5*(3.* ( pow((r+dx),2)*fields_vect[0][ind_space+1] - pow((r-dx),2)*fields_vect[0][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))
            -0.5*(3.* ( pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))
            );
}

double model1_charvar_PsiMinus(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    return(-fields_vect[0][ind_space]*fields_vect[1][ind_space]
            -Dx[0](fields_vect[0],ind_space,dx)/2.
            -Dx[0](fields_vect[1],ind_space,dx)/2.
            +0.5*(3.* ( pow((r+dx),2)*fields_vect[0][ind_space+1] - pow((r-dx),2)*fields_vect[0][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))
            -0.5*(3.* ( pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))
            );
}

//------------------- MODEL 1: CHARACTERISTIC VARIABLE, CHI RESCALING, HYPERBOLOIDAL FOLIATION -------------------//
//fields_vect [0]-> PsiPlus, [1]->PsiMinus, [2]->Psi

double model1_charvar_compactified_Psi_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    double s = param[0];
    double A = -r*r+s*s;
    double G = pow(r,4)+pow(s,4)+r*r*s*s*(s*s-2.);
    
    return ( 0.5*fields_vect[1][ind_space]+A*fields_vect[0][ind_space]/2/sqrt(G));
}

double model1_charvar_compactified_PsiPlus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    double s = param[0];
    
    double A = -r*r+s*s;
    double B = pow(r,7)-pow(r,5)*s*s*(2+s*s)+r*pow(s,6)*(2+s*s);
    double C = pow(r,3)+r*s*s*(-2+s*s);
    double D = pow(r,4) - 4.*r*r*s*s-pow(s,4);
    double G = pow(r,4)+pow(s,4)+r*r*s*s*(s*s-2.);
    double M = r*r+s*s;
    double N = 0.5*(-2.+M*s*s/G);
    double O = -2.*A*A+pow(s,8)/G;
    return (B*fields_vect[0][ind_space]/D/G
            +M*fields_vect[1][ind_space]*(C+A*s*s*fields_vect[0][ind_space])/sqrt(G)/D
            -M*pow(s,4)*Dx[0](fields_vect[1],ind_space,dx)/2./D/sqrt(G)
            +M*pow(s,4)*(3.* ( pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3))) / (2.*D)/sqrt(G)
            +O*Dx[0](fields_vect[0],ind_space,dx)/(2.*D)
            -pow(s,8)*(3.* ( pow((r+dx),2)*fields_vect[0][ind_space+1] - pow((r-dx),2)*fields_vect[0][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))/(2.*D*G)
            );
}

double model1_charvar_compactified_PsiMinus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double r = dmin+dx*double(ind_space-gl);
    double s = param[0];
    
    double G = pow(r,4)+pow(s,4)+r*r*s*s*(s*s-2.);
    double M = r*r+s*s;
    double N = 0.5*(-2.+M*s*s/G);
    
    return (M*fields_vect[1][ind_space]*(r-s*s*fields_vect[0][ind_space])/G
            +N*Dx[0](fields_vect[1],ind_space,dx)
            -M*s*s*(3.* ( pow((r+dx),2)*fields_vect[1][ind_space+1] - pow((r-dx),2)*fields_vect[1][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))/(2*G)
            -M*Dx[0](fields_vect[0],ind_space,dx)/(2*sqrt(G))
            +M*(3.* ( pow((r+dx),2)*fields_vect[0][ind_space+1] - pow((r-dx),2)*fields_vect[0][ind_space-1] )/(pow((r+dx),3)-pow((r-dx),3)))/(2*sqrt(G))
            );
}

// ----------- // MODEL 2 // ----------- //
double model2_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +param[0]*(pow(fields_vect[4][ind_space],2)-pow(fields_vect[3][ind_space],2)) );
}

double model2_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +param[0]*(pow(fields_vect[1][ind_space],2)-pow(fields_vect[0][ind_space],2)) );
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
    return (3.* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
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

// ----------- // DG model // ----------- //

// version 1
double DG_Model1_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3.* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            -2.*fields_vect[3][ind_space]/pow(1.+x*x,0.5)
            +sin(log(1.+pow(1.+x*x,0.5)*fields_vect[2][ind_space]))*((1.+x*x)*(-fields_vect[0][ind_space]+fields_vect[1][ind_space]+x*fields_vect[2][ind_space]))
            /pow(1.+x*x,1.5)/(1+pow(1.+x*x,0.5)*fields_vect[2][ind_space])
            //+ pow(fields_vect[3][ind_space],2)-pow(fields_vect[4][ind_space])
            );
}

// version 2
double DG_Model2_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3.* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            -2.*fields_vect[3][ind_space]/pow(1.+x*x,0.5)
            -sin(log(1.+pow(1.+x*x,0.5)*fields_vect[2][ind_space]))*((1.+x*x)*(-fields_vect[0][ind_space]+fields_vect[1][ind_space]+x*fields_vect[2][ind_space]))
            /pow(1.+x*x,1.5)/(1.+pow(1.+x*x,0.5)*fields_vect[2][ind_space])
            //+ pow(fields_vect[3][ind_space],2)-pow(fields_vect[4][ind_space],2)
            );
}
// ----------- // MODEL 3 Bis // ----------- //
// field [0]P1, [1]PHI1, [2]phi1, [3]P2, [4]PHI, [5]phi
double model3bis_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[2][ind_space] + param[0] * fields_vect[5][ind_space])/pow(param[0],2)*
            (pow(x*fields_vect[1][ind_space]+fields_vect[3][ind_space],2)+pow(x*fields_vect[4][ind_space]+fields_vect[5][ind_space],2)-pow(x*fields_vect[0][ind_space],2)-pow(x*fields_vect[3][ind_space],2)+2*(fields_vect[2][ind_space]+fields_vect[5][ind_space]))
            );
}

double model3bis_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (3* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) + 
            (fields_vect[5][ind_space] - param[0] * fields_vect[2][ind_space])/pow(param[0],2)*
            (pow(x*fields_vect[1][ind_space]+fields_vect[3][ind_space],2)+pow(x*fields_vect[4][ind_space]+fields_vect[5][ind_space],2)-pow(x*fields_vect[0][ind_space],2)-pow(x*fields_vect[3][ind_space],2)+2*(fields_vect[2][ind_space]+fields_vect[5][ind_space]))
            );
}


// ----------- // MODEL 3 Bis General version // ----------- //

// ----K=1+R^2---- //

// field [0]P1, [1]PHI1, [2]phi1, [3]P2, [4]PHI, [5]phi
double model3bisGen_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double A = param[0];
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[1], ind_space, dx) - (1/(A*A)) *2*(1 + R*R)*(2*R*fields_vect[1][ind_space] + 
    fields_vect[2][ind_space])*(pow(fields_vect[2][ind_space],2) + 
    pow(fields_vect[1][ind_space],2)) - (1/(A*A))*(fields_vect[2][ind_space] + A*fields_vect[5][ind_space])*(-pow((1 + R*R),2)*pow(fields_vect[0][ind_space],2) + pow((1 + R*R)*fields_vect[1][ind_space] + 2*R*fields_vect[2][ind_space],2) - pow((1 + R*R),2)*pow(fields_vect[3][ind_space],2) + pow((1 + R*R)*fields_vect[4][ind_space] + 2*R*fields_vect[5][ind_space],2))
        );
}

double model3bisGen_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double A = param[0];
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[4], ind_space, dx) - (1./(A*A))*2.* (1. + R*R) *(2.* R *fields_vect[4][ind_space] + 
    fields_vect[5][ind_space]) *(pow(fields_vect[2][ind_space],2) + 
    pow(fields_vect[5][ind_space],2)) - (1./(A*A))*(-A*fields_vect[2][ind_space] + 
    fields_vect[5][ind_space])*(-pow((1. + R*R),2) *pow(fields_vect[0][ind_space],2) + 
    pow((1.+R*R)*fields_vect[1][ind_space]+2.*R*fields_vect[2][ind_space],2) - pow((1. + R*R),2)*pow(fields_vect[3][ind_space],2) + 
    pow((1. + R*R)*fields_vect[4][ind_space] + 2.*R*fields_vect[5][ind_space],2)));
}

// EoM for the 'rescaled' fields
// field [0]KP1, [1]KPHI1, [2]Kphi1, [3]KP2, [4]KPHI, [5]Kphi

double model3bisGen_KPI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[1],ind_space,dx)
        -(fields_vect[2][ind_space] + param[0] * fields_vect[5][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2)));
}

double model3bisGen_KPI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[4],ind_space,dx)
        -(-fields_vect[5][ind_space] + param[0] * fields_vect[2][ind_space])/pow(param[0],2)*
            (pow(fields_vect[1][ind_space],2)+pow(fields_vect[4][ind_space],2)-pow(fields_vect[0][ind_space],2)-pow(fields_vect[3][ind_space],2)));
}

// ----K=(1+R^2)^0.5---- //

double model3bisGen_PI1_sqr(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double A = param[0];
    double K = pow(R*R+1.,0.5);
    double KPrime = R/pow(R*R+1.,0.5);
    double KDoublePrime = 1./pow(R*R+1.,1.5);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[1],ind_space,dx) -
            1./(A*A)*K*(2*KPrime*fields_vect[1][ind_space] + KDoublePrime*fields_vect[2][ind_space]) 
            *(pow(fields_vect[2][ind_space],2)+pow(fields_vect[5][ind_space],2)) 
            - 1./(A*A)*(fields_vect[2][ind_space]+A*fields_vect[5][ind_space]) 
            *(-K*K* pow(fields_vect[0][ind_space],2) + pow((K*fields_vect[1][ind_space]+KPrime*fields_vect[2][ind_space]),2) - 
    K*K*pow(fields_vect[3][ind_space],2) + pow((K*fields_vect[4][ind_space] + KPrime*fields_vect[5][ind_space]),2))
);
}


double model3bisGen_PI2_sqr(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double A = param[0];
    double K = pow(R*R+1.,0.5);
    double KPrime = R/pow(R*R+1.,0.5);
    double KDoublePrime = 1./pow(R*R+1.,1.5);
    //cout<<"x in eq for PI "<<x<<endl;
    //cout<< " dt dx gl ord espilon "<<dt<<" "<<dx<<" "<<gl<<" "<<ord<<" "<<epsilon<<endl;
    return (-Dx[0](fields_vect[4], ind_space, dx) - 
            (1/(A*A))*K *(2* KPrime* fields_vect[4][ind_space] + KDoublePrime *fields_vect[5][ind_space])* 
            (pow(fields_vect[2][ind_space],2) + pow(fields_vect[5][ind_space],2)) - 
            1/(A*A)*(-A*fields_vect[2][ind_space] + fields_vect[5][ind_space])*(-K*K* pow(fields_vect[0][ind_space],2) 
            + pow((K*fields_vect[1][ind_space] +KPrime*fields_vect[2][ind_space]),2) - 
            K*K*pow(fields_vect[3][ind_space],2) + pow((K*fields_vect[4][ind_space]+KPrime*fields_vect[5][ind_space]),2)));
}


// MODEL 3 HYBRID //
// we mix the model 3 equations with the dumped model 3bis equation
// we take the version of model 3bis with K=(R^2+1)^(0.5)
// dumping coefficient in front of the model3bis equation: R^2/(R^2+1)^0.5

double model3Hyb_PI1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double x = R;
    double A = param[0];
    double K = R+1.;
    double KPrime = 1.;
    double KDoublePrime = 0.;
    return (//wave equation part
            (1./(R+1.))*
            (3.*(pow((x+dx),2)*fields_vect[1][ind_space+1]-pow((x-dx),2)*fields_vect[1][ind_space-1])/(pow((x+dx),3)-pow((x-dx),3)))
            
            // model 3bis part
            
            +pow(R,2)/(pow(R,2)+1.)*
            (-Dx[0](fields_vect[1],ind_space,dx) -
            (2. *KPrime *fields_vect[1][ind_space] + KDoublePrime *fields_vect[2][ind_space])/K
            - 1./(A*A*K)*(fields_vect[2][ind_space]+A*(fields_vect[5][ind_space])) 
            *(-K*K* pow(fields_vect[0][ind_space],2) + pow((K*fields_vect[1][ind_space]+KPrime*fields_vect[2][ind_space]),2) - 
            K*K*pow(fields_vect[3][ind_space],2) + pow(K*fields_vect[4][ind_space] + KPrime*fields_vect[5][ind_space],2)))
            );
}

double model3Hyb_PI2(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double R = dmin+dx*(ind_space-gl);
    double x = R;
    double A = param[0];
    double K = R+1.;
    double KPrime = 1.;
    double KDoublePrime = 0.;
    return (// wave equation part
            (1./(R+1.))*
            (3.*(pow((x+dx),2)*fields_vect[4][ind_space+1]-pow((x-dx),2)*fields_vect[4][ind_space-1])/(pow((x+dx),3)-pow((x-dx),3)))
            // model 3bis part
            
            +pow(R,2)/(pow(R,2)+1.)*
            (-Dx[0](fields_vect[4], ind_space, dx) - 
            (2.* KPrime* fields_vect[4][ind_space] + KDoublePrime* fields_vect[5][ind_space])/K 
            -1./(A*A*K)*(-A*fields_vect[2][ind_space] + fields_vect[5][ind_space])*(-K*K* pow(fields_vect[0][ind_space],2) 
            + pow((K*fields_vect[1][ind_space] +KPrime*fields_vect[2][ind_space]),2) - 
            K*K*pow(fields_vect[3][ind_space],2) + pow((K*fields_vect[4][ind_space]+KPrime*fields_vect[5][ind_space]),2)))
            );
}



//------------------- MODEL 3: CHARACTERISTIC VARIABLE, CHI RESCALING, HYPERBOLOIDAL FOLIATION -------------------//

double model3_charvar_compactified_PsiOne_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double s = param[0];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2.+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    return ( 0.5*fields_vect[1][ind_space] + N*fields_vect[0][ind_space]/2./sqrt(J) );
}

double model3_charvar_compactified_PsiTwo_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double s = param[0];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2.+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    return ( 0.5*fields_vect[4][ind_space] + N*fields_vect[3][ind_space]/2./sqrt(J) );
}

double model3_charvar_compactified_PsiOnePlus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double x = r;
    double s = param[0];
    double A = param[3];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    double P = (pow(r,7)-pow(r,5)*s*s*(2.+s*s)+r*pow(s,6)*(2.+s*s));

    return ( P*fields_vect[0][ind_space]/(J*K)
            -O*s*s*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(sqrt(J)*K)
            +(fields_vect[2][ind_space]/A+fields_vect[5][ind_space])*L*L*M*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*J*K)
            +sqrt(J)*M/(2.*K)*(-Dx[0](fields_vect[1],ind_space,dx) + 3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)))
            +(-2.*J*N*N+pow(s,8))*Dx[0](fields_vect[0],ind_space,dx)/(2*J*K)
            -pow(s,8)* 3* ( pow((x+dx),2)*fields_vect[0][ind_space+1] - pow((x-dx),2)*fields_vect[0][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))/(2*J*K)
    );
}

double model3_charvar_compactified_PsiTwoPlus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double x = r;
    double s = param[0];
    double A = param[4];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    double P = (pow(r,7)-pow(r,5)*s*s*(2.+s*s)+r*pow(s,6)*(2.+s*s));
    return ( P*fields_vect[3][ind_space]/(J*K)
            -O*s*s*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*sqrt(J)*K)
            +(-fields_vect[2][ind_space]+fields_vect[5][ind_space]/A)*L*L*M*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*J*K)
            +sqrt(J)*M/(2.*K)*(-Dx[0](fields_vect[4],ind_space,dx) + 3.* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)))
            +(-2.*J*N*N+pow(s,8))*Dx[0](fields_vect[3],ind_space,dx)/(2*J*K)
            -pow(s,8)* 3.* ( pow((x+dx),2)*fields_vect[3][ind_space+1] - pow((x-dx),2)*fields_vect[3][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))/(2*J*K)
    );
}

double model3_charvar_compactified_PsiOneMinus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double x = r;
    double s = param[0];
    double A = param[3];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    return ( -M*s*s*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/J
            +O*s*s*(fields_vect[2][ind_space]/A+fields_vect[5][ind_space])*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*pow(J,3./2.))
            -(2*J+O)*Dx[0](fields_vect[1],ind_space,dx)/(2.*J)
            +O/(2*J)*3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))
            +M/(2*sqrt(J))*(-Dx[0](fields_vect[0],ind_space,dx)+3* ( pow((x+dx),2)*fields_vect[0][ind_space+1] - pow((x-dx),2)*fields_vect[0][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)))
    );
}

double model3_charvar_compactified_PsiTwoMinus_Chi(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    //fields_vect [0]-> Psi1Plus, [1]->Psi1Minus, [2]->Psi1, [3]-> Psi2Plus, [4]->Psi2Minus, [5]->Psi2
    double r = dmin+dx*double(ind_space-gl);
    double x = r;
    double s = param[0];
    double A = param[4];
    double K = pow(r,4)-4*r*r*s*s-pow(s,4);
    double J = pow(r,4)+pow(s,4)+r*r*s*s*(-2+s*s);
    double L = -r*r*s+pow(s,3);
    double M = r*r+s*s;
    double N = -r*r+s*s;
    double O = pow(r,4)-pow(s,4);
    return ( -M*s*s*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*J)
            -O*s*s*(fields_vect[2][ind_space]-fields_vect[5][ind_space]/A)*(fields_vect[1][ind_space]*fields_vect[0][ind_space]+fields_vect[4][ind_space]*fields_vect[3][ind_space])/(A*pow(J,3./2.))
            -(2*J+O)*Dx[0](fields_vect[4],ind_space,dx)/(2.*J)
            +O/(2.*J)*3.* ( pow((x+dx),2)*fields_vect[4][ind_space+1] - pow((x-dx),2)*fields_vect[4][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))
            +M/(2.*sqrt(J))*(-Dx[0](fields_vect[3],ind_space,dx)+3* ( pow((x+dx),2)*fields_vect[3][ind_space+1] - pow((x-dx),2)*fields_vect[3][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)))
    );
}
// Field [0] Pi,[1]Phi,[2]Psi,[3]F 
//------------------- Curved model for the wave equation -------------------//

// H=0, G=1/F
double curved_wave_eq_spherical_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    
    return ( pow(fields_vect[3][ind_space],2)*(3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))) );
}

// H=0, G=F/(2-F)
double curved_wave_eq_spherical_PI_1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    
    return ( -(-2.+pow(fields_vect[3][ind_space],2))*(3.* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3))) );
}

double curved_wave_eq_spherical_F(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    
    return(-d*(fields_vect[3][ind_space]-1)-K*pow(  fields_vect[0][ind_space]+fields_vect[3][ind_space]*fields_vect[2][ind_space],2)-fields_vect[3][ind_space]*Dx[0](fields_vect[3],ind_space,dx));
}

double curved_wave_eq_spherical_F1(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    if(fields_vect[3][ind_space]>0)
    {
        return(-d*(fields_vect[3][ind_space]-1)-K*pow(  fields_vect[0][ind_space]+fields_vect[3][ind_space]*fields_vect[2][ind_space],2)-fields_vect[3][ind_space]*Dx[0](fields_vect[3],ind_space,dx));
    }
    if(fields_vect[3][ind_space]<=0)
    {
        return(-d/(pow(fields_vect[3][ind_space],3)-1)*fields_vect[3][ind_space]-K*pow(  fields_vect[0][ind_space]+fields_vect[3][ind_space]*fields_vect[2][ind_space],2)-fields_vect[3][ind_space]*Dx[0](fields_vect[3],ind_space,dx));
    }
    return 0;
}

double curved_wave_eq_spherical_F1_parity_condition(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    if(fields_vect[3][ind_space]>0)
    {
        return(-d*(fields_vect[3][ind_space]-1)-K*pow(  fields_vect[0][ind_space]+fields_vect[3][ind_space]*fields_vect[2][ind_space]*fields_vect[2][ind_space],2)-fields_vect[3][ind_space]*Dx[0](fields_vect[3],ind_space,dx)*Dx[0](fields_vect[3],ind_space,dx));
    }
    if(fields_vect[3][ind_space]<=0)
    {
        return(-d/(pow(fields_vect[3][ind_space],3)-1)*fields_vect[3][ind_space]-K*pow(  fields_vect[0][ind_space]+fields_vect[3][ind_space]*fields_vect[2][ind_space]*fields_vect[2][ind_space],2)-fields_vect[3][ind_space]*Dx[0](fields_vect[3],ind_space,dx)*Dx[0](fields_vect[3],ind_space,dx));
    }
    return 0;
}

double curved_wave_eq_spherical_F_simple(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    if(fields_vect[3][ind_space]>0)
    {
        return(-d*(fields_vect[3][ind_space]-1)*fields_vect[3][ind_space]-K*(fields_vect[2][ind_space]));
    }
    if(fields_vect[3][ind_space]<=0)
    {
        return(-d/(pow(fields_vect[3][ind_space],3)-1)*fields_vect[3][ind_space]-K*(fields_vect[2][ind_space]));
    }
    return 0;
}

double curved_wave_eq_spherical_F_simple_der(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    if(fields_vect[3][ind_space]>0)
    {
        return(-d*(fields_vect[3][ind_space]-1)*fields_vect[3][ind_space]-K*(pow(fields_vect[0][ind_space],2)+pow(fields_vect[1][ind_space],2)));
    }
    if(fields_vect[3][ind_space]<=0)
    {
        return(-d/(pow(fields_vect[3][ind_space],3)-1)*fields_vect[3][ind_space]-K*(pow(fields_vect[0][ind_space],2)+pow(fields_vect[1][ind_space],2)));
    }
    return 0;
}

//------------------- Curved model for the model 1 equations -------------------//

double curved_m1_spherical_PI(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double x = dmin+dx*(ind_space-gl);
    
    return ( -pow(fields_vect[0][ind_space],2)+pow(fields_vect[3][ind_space],2)*(pow(fields_vect[1][ind_space],2)+(3* ( pow((x+dx),2)*fields_vect[1][ind_space+1] - pow((x-dx),2)*fields_vect[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)))) );
}

double curved_m1_spherical_F_simple(int ind_field,int ind_space,std::vector<std::vector<double>> fields_vect,double dx,double dmin,std::vector<double> param, double t,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,double dt, int gl)
{
    double K = param[0];
    double d = param[3];
    double x = dmin+dx*(ind_space-gl);
    if(fields_vect[3][ind_space]>0)
    {
        return(-d*(fields_vect[3][ind_space]-1)*fields_vect[3][ind_space]-K*(fields_vect[2][ind_space]));
    }
    if(fields_vect[3][ind_space]<=0)
    {
        return(-d/(pow(fields_vect[3][ind_space],3)-1)*fields_vect[3][ind_space]-K*(fields_vect[2][ind_space]));
    }
    return 0;
}

// ---------- INITIAL DATA AND INITIALIZATION OF VECTORS --------------- //

double initial_line(double x,double init_param, vector<double> parameter) 
{
    return(x);
}

double initial_unity(double x,double init_param, vector<double> parameter) 
{
    return(1.);
}

double initial_unity_K(double x,double init_param, vector<double> parameter) 
{
    return(1.*(1.+x*x));
}

double initial_null(double x,double init_param, vector<double> parameter) 
{
    double a = 0.;
    return(a);
}


double initial_gauss(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    return( init_param*exp(-pow(((x-mu)*dev_std),2)) );
}

double initial_gauss_m1(double x,double init_param, vector<double> parameter) 
{
    double A = parameter[0];
    double dev_std = parameter[1];
    double mu = parameter[2];
    return( log(1+A*init_param*exp(-pow(((x-mu)*dev_std),2)) ));
}

double initial_gauss_m2(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    double second_init_param = parameter[3];
    return( second_init_param*exp(-pow(((x-mu)*dev_std),2)) );
}

double initial_gauss_m2_PHI(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    double second_init_param = parameter[3];
    return( second_init_param*(-2)*(x-mu)*dev_std*dev_std*exp(-pow(((x-mu)*dev_std),2)) );
}


double initial_gauss_K(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    double k = 1.+x*x;
    return( init_param*exp(-pow(((x-mu)*dev_std),2))*(1+x) );
}

double initial_gauss_m3(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    return(sqrt( 1-pow(init_param*exp(-pow(((x-mu)*dev_std),2)),2)) );  
}

double initial_gauss_m3_PHI(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    double second_init_param = parameter[3];
    return(    0.5*pow(( 1-pow(init_param*exp(-pow(((x-mu)*dev_std),2)),2)),-0.5)*(-2.)*init_param*exp(-pow(((x-mu)*dev_std),2))*(-2)*(x-mu)*dev_std*dev_std);
}

double initial_gauss_PHI(double x,double init_param, vector<double> parameter) 
{
    double dev_std = parameter[1];
    double mu = parameter[2];
    return( -2.*(x-mu)*dev_std*dev_std*init_param*exp(-pow(((x-mu)*dev_std),2)) );
}

double initial_gauss_PHI_compactified(double x,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
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

double initial_gauss_phi_compactified(double x,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
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


// ------ hyperboloidal compactification and Chi function rescaling  for the wave equation------ //

double initial_gauss_PI_compactified_Chi1_we(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
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
        return( (-A*exp(-C)+A*exp(B)+2.*A*C*exp(-C)-2.*A*D*ds_sq*exp(B))*Chi/2./R);
    }
    else if (r==0)
    {
        return(0);
    }
    else
    {
        return(1./2.*(-A*exp(-C)+2.*A*C*exp(-C)));
    }
}


double initial_gauss_PHI_compactified_Chi_we(double x,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double rate_of_square = pow(x,2)/pow(s,2);
    double R = x/(1.-rate_of_square);
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

double initial_gauss_PI_compactified_Chi_we_null_time_derivative(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = -dev_std*dev_std*(-r+2.*R)*(-r+2.*R);
    double C = dev_std*dev_std*r*r;
    double G = (-r+2.*R);
    if(r!=0 && r!=s)
    {
        return( A*Chi*(exp(B) - exp(-C) + 
       2.*pow(dev_std,2)*
        (-pow(G,2)*exp(B) + exp(4*pow(dev_std,2)*R*(-r + R)-dev_std*dev_std*(-r+2.*R)*(-r+2.*R))*pow(r,2)))/
   (2.*R));
    }
    else if (r==0)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PHI_compactified_Chi_we_null_time_derivative(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = exp(-pow(dev_std,2)*pow(r-2.*R,2));
    double C = r-2.*R;
    double D = -r+R;
    double G = dev_std*dev_std*r*r;
    double H = exp(-dev_std*dev_std*(r*r+pow(r-2.*R,2)));
    
    if(r!=s && r!=0)
    {        
        return (A*H*(-exp(G)*(-Chi*r + C*ChiPrime*R + 
            Chi*R*(1. + 2.*G + 8.*D*pow(dev_std,2)*R))) + 
       (exp(pow(C,2)*pow(dev_std,2)-dev_std*dev_std*(r*r+pow(r-2.*R,2)))*(ChiPrime*r*R + Chi*(D - 2.*G*R)))/
   (2.*pow(R,2)) );
    }    
    else if(r==0)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_phi_compactified_Chi_we_null_time_derivative(double r,double init_param, vector<double> parameter) 
{
    double A = init_param;
    double s = parameter[0];
    double dev_std = parameter[1];
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);

    if(r!=s & r!= 0)
    {
        return( init_param*Chi/(2.*R)*  (  exp(-pow(dev_std*r,2))*r-exp(-pow(dev_std*(r-2.*R),2))*(r-2.*R))  );
    }
    else if (r==s)
    {
        return(0.5*A*exp(-pow(dev_std*r,2))*s);
    }
    else
    {
        return(A);
    }
}

double initial_gauss_phi_compactified_Chi_we_null_field(double r,double init_param, vector<double> parameter) 
{
    double A = init_param;
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);

    if(r!=s & r!=0)
    {
        return( A*(-exp(-ds*ds*r*r)+exp(-ds*ds*pow(r-2.*R,2)))*Chi/R);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PHI_compactified_Chi_we_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = exp(-pow(dev_std,2)*pow(r-2.*R,2));
    double C = r-2.*R;
    double D = -r+R;
    double G = dev_std*dev_std*r*r;
    double H = exp(-dev_std*dev_std*(r*r+pow(r-2.*R,2)));
    
    if(r!=s && r!=0)
    {        
        return (A*(pow(exp(1.),pow(C,2)*pow(ds,2)-dev_std*dev_std*(r*r+pow(r-2.*R,2)))*
        (Chi - ChiPrime*R + 2.*Chi*pow(ds,2)*r*R) + H*
       pow(exp(1.),G)*(ChiPrime*R + Chi*(-1. + 2.*C*pow(ds,2)*R)))/pow(R,2));
    }    
    else if(r==0)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PI_compactified_Chi_we_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = -dev_std*dev_std*(-r+2.*R)*(-r+2.*R);
    double C = dev_std*dev_std*r*r;
    double G = (-r+2.*R);
    if(r!=0 && r!=s)
    {
        return((2*A*Chi*pow(ds,2)*(-(pow(exp(1.),B)*G) - r/pow(exp(1.),C)))/R);
    }
    else if (r==0)
    {
        return(-4.*A*ds*ds);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_phi_compactified_Chi1_we(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
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

// -------- MODEL1 char variable, flat space -------- //

double initial_gauss_PsiPlus_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double a = init_param;
    double b = 1;
    double ds = 1;
    return(-(2*a*b*ds*ds*r)/(a*b+exp(ds*ds*r*r)));
}

double initial_gauss_PsiMinus_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double a = init_param;
    double b = 1;
    double ds = 1;
    return((2*a*b*ds*ds*r)/(a*b+exp(ds*ds*r*r)));
}

double initial_gauss_Psi_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double a = init_param;
    double b = 1;
    double ds = 1;
    return(log(1+a*b*exp(-ds*ds*r*r)));
}

// ------ hyperboloidal compactification and Chi function rescaling for the MODEL 1 equations  with char variables------ //

double initial_gauss_Psi_compactified_Chi_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    if( r!=0 && r!=s)
    {
        return(log(1.+A*B*(-exp(-ds_sq*r*r)+exp(-ds_sq*pow(r-2.*R,2)))/R)*Chi);
    }
    
    else if (r==s)
    {
        return(B*A/2.*exp(-ds_sq*r*r)*r);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PsiPlus_compactified_Chi_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/((double)1.-rate_of_square);
    double Chi = pow((double)1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    double K = exp(-ds_sq*pow(r-2.*R,2));
    double J = exp( ds_sq*r*r);
    //cout<<"K: "<<K<<endl;
    if( r!=0 && r!=s)
    {
        return(A*B*(-r+J*K*(r-(double)4.*ds_sq*pow(r-(double)2.*R,2)*R)) *Chi*Chi
                /(R*(A*B*(r-J*K*(r-(double)2.*R))+(double)2.*J*R)) );
    }
    else if (r==s)
    {
        return(-B*A/(double)2.*r*exp(-ds_sq*r*r));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PsiMinus_compactified_Chi_charvar_m1(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    double K = exp( -pow(dev_std*(r - 2. *R),2));
    double J = exp( pow(dev_std*r,2));
    if( r!=0 && r!=s)
    {
        return(A*B*(-J*K*(r-2.*R)+(r+(-2.+4.*ds_sq*r*r)*R))*Chi
                /(R*(A*B*(r-J*K*(r-2.*R))+2.*J*R))
        );
    }
    else if (r==s)
    {
        return(B*A/2.*exp(-ds_sq*r*r)*(-2.+4.*ds_sq*r*r));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi_compactified_Chi_charvar_m1_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    if( r!=0 && r!=s)
    {
        return(log(1.+A*B*(-exp(-ds_sq*r*r)+exp(-ds_sq*pow(r-2.*R,2)))/R)*Chi);
    }
    
    else if (r==s)
    {
        return(-((A*B*s)/(pow(exp(1),pow(ds,2)*pow(s,2))*sqrt(pow(s,2)))));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_PsiPlus_compactified_Chi_charvar_m1_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
        double ds = dev_std;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/((double)1.-rate_of_square);
    double Chi = pow((double)1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    double num;
    double den;
    //cout<<"K: "<<K<<endl;
    if( r!=0 && r!=s)
    {
        num = A*B*pow(Chi,2)*(1 + pow(exp(1),4*pow(ds,2)*(r - R)*R)*(-1 + 4*pow(ds,2)*(r - 2*R)*R));
        den = R*(A*B*(-1 + pow(exp(1),4*pow(ds,2)*(r - R)*R)) + pow(exp(1),pow(ds,2)*pow(r,2))*R);
        return(num/den);
    }
    else if (r==s)
    {
        return((A*B*(2*pow(s,4) + pow(s,4)*(-2 + pow(s,2))))/(pow(exp(1),pow(ds,2)*pow(s,2))*pow(s,6)));
    }
    else
    {
        return(-4.* A* B* ds*ds);
    }
}

double initial_gauss_PsiMinus_compactified_Chi_charvar_m1_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = 1.;
    double num;
    double den;
    
    if( r!=0 && r!=s)
    {
        num = A*B*Chi*(-1 + pow(exp(1),4*pow(ds,2)*(r - R)*R) - 4*pow(ds,2)*r*R);
        den = R*(A*B*(-1 + pow(exp(1),4*pow(ds,2)*(r - R)*R)) + pow(exp(1),pow(ds,2)*pow(r,2))*R); 
        return(num/den);
    }
    else if (r==s)
    {
        return((A*B*sqrt(pow(s,2))*(pow(s,2) + pow(s,2)*(-1 + 4*pow(ds,2)*pow(s,2))))/
   (pow(exp(1),pow(ds,2)*pow(s,2))*pow(s,4)));
    }
    else                        
    {
        return(-4.* A* B* ds*ds);
    }
}

// ------ hyperboloidal compactification and Chi function rescaling for the MODEL 3 equations  with char variables------ //

// initial data for F(x)=x*Aexp(-x^2)

double initial_gauss_Psi1_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = 1.;
    if( r!=0 && r!=s)
    {
        return(C*sin(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C)*Chi);
    }
    
    else if (r==s)
    {
        return(A/2.*exp(-ds_sq*r*r)*r);
    }
    else
    {
        return(C*sin(log(1.+A)/C));
    }
}

double initial_gauss_Psi2_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = parameter[5];
    if( r!=0 && r!=s)
    {
        return(C*(-1+cos(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C))*Chi);
    }
    
    else if (r==s)
    {
        return(0);
    }
    else
    {
        return(C*(-1+cos(log(1.+A)/C)));
    }
}

double initial_gauss_Psi1Plus_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = 1.;
    double K = exp( -pow(dev_std*(r - 2. *R),2));
    double J = exp( pow(dev_std*r,2));
    double L = r-2.*R;
    if( r!=0 && r!=s)
    {
        return(-A*cos(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C)
            *(r+K*J*(-r+4.*ds_sq*L*L*R))*Chi*Chi 
            /(R*(A*r+2.*J*R)+R*K*(-A*J*r+2.*A*J*R))
        );
    }
    
    else if (r==s)
    {
        return(-0.5*A*exp(-ds_sq*r*r)*r);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi2Plus_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = 1.;
    double K = exp( -pow(dev_std*(r - 2. *R),2));
    double J = exp( pow(dev_std*r,2));
    double L = r-2.*R;
    if( r!=0 && r!=s)
    {
        return(A*sin(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C)
            *(r+K*J*(-r+4.*ds_sq*L*L*R))*Chi*Chi 
            /(R*(A*r+2.*J*R)+R*K*(-A*J*r+2.*A*J*R))
        );
    }
    
    else if (r==s)
    {
        return(-0.5*A*exp(-ds_sq*r*r)*r);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi1Minus_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = 1.;
    double K = exp( -pow(dev_std*(r - 2. *R),2));
    double J = exp( pow(dev_std*r,2));
    double L = r-2.*R;
    if( r!=0 && r!=s)
    {
        return(A*cos(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C)
            *(-J*L*K+(r+(-2.+4.*ds_sq*r*r)*R))*Chi 
            /(R*(A*r+2.*J*R)+R*K*(-A*J*r+2.*A*J*R))
        );
    }
    
    else if (r==s)
    {
        return(A*exp(-ds_sq*r*r)*(-1+2*ds_sq*r*r));
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi2Minus_compactified_Chi_charvar_m3(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double C = 1.;
    double K = exp( -pow(dev_std*(r - 2. *R),2));
    double J = exp( pow(dev_std*r,2));
    double L = r-2.*R;
    if( r!=0 && r!=s)
    {
        return(A*sin(log(1.+(A*exp(-ds_sq*r*r)*r-A*exp(-ds_sq*pow(r-2.*R,2))*(r-2.*R))/(2.*R))/C)
            *(J*L*K-(r+(-2.+4.*ds_sq*r*r)*R))*Chi 
            /(R*(A*r+2.*J*R)+R*K*(-A*J*r+2.*A*J*R))
        );
    }
    
    else if (r==s)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

// --- F(x)=A*exp(-x^2)

double initial_gauss_Psi1_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1.-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return(C1*Chi*sin(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/
        R)/C1));
    }
    
    else if (r==s)
    {
        return(A *exp(-ds_sq*s*s));
    }
    else
    {
        return(0);
    }
}


double initial_gauss_Psi2_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double C2 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return(C2*Chi*(-1+cos(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/R)/C2)));
    }
    
    else if (r==s)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi1Plus_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double C2 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return((A*pow(Chi,2)*(1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))*
        (-1 + 4*pow(ds,2)*(r - 2*R)*R))*
     cos(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/
          R)/C1))/(R*(A*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))) + 
       pow(exp(1),pow(ds,2)*(pow(r,2) + pow(r - 2*R,2)) - pow(ds,2)*pow(r - 2*R,2))*R)));
    }
    
    else if (r==s)
    {
        return(A/pow(exp(1),pow(ds,2)*pow(s,2)));
    }
    else
    {
        return(-4.* A* ds_sq);
    }
}

double initial_gauss_Psi2Plus_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double C2 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return((A*pow(Chi,2)*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))*
        (1 + 4*pow(ds,2)*R*(-r + 2*R)))*
     sin(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/
          R)/C2))/(R*(A*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))) + 
       pow(exp(1),pow(ds,2)*pow(r,2))*R)));
    }
    
    else if (r==s)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

double initial_gauss_Psi1Minus_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return((A*Chi*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2)) - 4*pow(ds,2)*r*R)*
     cos(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/
          R)/C1))/(R*(A*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))) + 
       pow(exp(1),pow(ds,2)*pow(r,2))*R)));
    }
    
    else if (r==s)
    {
        return(4.* A *ds_sq* exp(-ds_sq *s*s)*s);
    }
    else
    {
        return(-4.* A* ds_sq);
    }
}

double initial_gauss_Psi2Minus_compactified_Chi_charvar_m3_null_field(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double ds = dev_std;
    double C1 = 1.;
    double C2 = 1.;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double Chi = pow(1.+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    if( r!=0 && r!=s)
    {
        return((A*Chi*(1 - pow(exp(1),pow(ds,2)*pow(r,2 - pow(ds,2)*pow(r - 2*R,2))) + 4*pow(ds,2)*r*R)*
     sin(log(1 + (A*(-pow(exp(1),-(pow(ds,2)*pow(r,2))) + pow(exp(1),-(pow(ds,2)*pow(r - 2*R,2)))))/
          R)/C2))/(R*(A*(-1 + pow(exp(1),pow(ds,2)*pow(r,2) - pow(ds,2)*pow(r - 2*R,2))) + 
       pow(exp(1),pow(ds,2)*pow(r,2))*R)));
    }
    
    else if (r==s)
    {
        return(0);
    }
    else
    {
        return(0);
    }
}

// ------ hyperboloidal compactification and Chi function rescaling for the MODEL 4 equations  with char variables------ //
double initial_gauss_Psi1_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double mu = parameter[2];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double A = init_param;
        return(A*exp(-ds_sq*pow(r-mu,2)));
}

double initial_gauss_Psi2_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double mu = parameter[2];
    double dev_std = parameter[1];
    double ds_sq = dev_std*dev_std;
    double A = init_param;
    double B = parameter[5];
        return(B*exp(-ds_sq*pow(r-mu,2)));
}
// if the PsiMinus' are null we can chose:
double initial_null_Psi1Minus_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double mu = parameter[2];
    double A = init_param;
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double RPrime = (1 + rate_of_square)/pow(1-rate_of_square,2);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    if(r!=s)
    {
        return(2.*A*exp(-ds_sq*pow(r-mu,2))*(2.*ds_sq*r*Chi-2.*ds_sq*mu*Chi+RPrime*ChiPrime)/(-1+2.*RPrime+Chi));
    }
    
    else 
    {
        return(-A*exp(-ds_sq*pow(s-mu,2)));
    }
}

double initial_null_Psi2Minus_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double mu = parameter[2];
    double A = init_param;
    double B = parameter[5];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double RPrime = (1 + rate_of_square)/pow(1-rate_of_square,2);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    if(r!=s)
    {
        return(4.*B*ds_sq*exp(-ds_sq*pow(r-mu,2))*(r-mu)*pow(Chi,2)/(-1.+2.*RPrime+Chi));
    }
    
    else 
    {
        return(B*ds_sq*exp(-ds_sq*pow(r-mu,2))*pow(s,2)*(s-mu));
    }
    
}

double initial_gauss_Psi1Plus_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double mu = parameter[2];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double RPrime = (1 + rate_of_square)/pow(1-rate_of_square,2);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    
    if(r!=s)
    {
        return(-2.*A*exp(-ds_sq*pow(r-mu,2))*(2.*ds_sq*r*Chi-2.*ds_sq*mu*Chi+RPrime*ChiPrime)/(-1+2.*RPrime+Chi));
    }
    
    else 
    {
        return(A*exp(-ds_sq*pow(s-mu,2)));
    }
}

double initial_gauss_Psi2Plus_compactified_Chi_charvar_m4(double r,double init_param, vector<double> parameter) 
{
    double s = parameter[0];
    double dev_std = parameter[1];
    double mu = parameter[2];
    double ds_sq = dev_std*dev_std;
    double rate_of_square = pow(r,2)/pow(s,2);
    double R = r/(1-rate_of_square);
    double RPrime = (1 + rate_of_square)/pow(1-rate_of_square,2);
    double Chi = pow(1+R*R,0.5);
    double ChiPrime = R /Chi;
    double A = init_param;
    double B = parameter[5];
    double C = 1.;
    
    if(r!=s)
    {
        return(-4.*B*ds_sq*exp(-ds_sq*pow(r-mu,2))*(r-mu)*pow(Chi,2)/(-1.+2.*RPrime+Chi));
    }
    
    else 
    {
        return(-B*ds_sq*exp(-ds_sq*pow(r-mu,2))*pow(s,2)*(s-mu));
    }
}
// ------  ------  ------  ------  ------  ------  ------  ------  ------  ------

// Model 3bis generalized initial functions with K=1+x^2 //

double initial_gauss_Psi1_m3bis_gen(double R,double init_param, vector<double> parameter) 
{
    double C1= parameter[0];
    double ds = parameter[1];
    double mu = parameter[2];
    double ds_sq = ds*ds;
    double A = init_param;
    double K = R*R+1;
    return((C1*sin(log(1. + (A*R)/exp(R*R*(ds_sq))/C1))/K));
}

double initial_gauss_Psi2_m3bis_gen(double R,double init_param, vector<double> parameter) 
{
    double C1= parameter[0];
    double ds = parameter[1];
    double mu = parameter[2];
    double ds_sq = ds*ds;
    double A = init_param;
    double K = R*R+1;
    
    return((C1*cos(log(1. + (A*R)/exp(R*R*(ds_sq))/C1))/K));
}

double initial_gauss_PHI1_m3bis_gen(double R,double init_param, vector<double> parameter) 
{
    double R_sq = R*R;
    double C1= parameter[0];
    double ds = parameter[1];
    double mu = parameter[2];
    double ds_sq = ds*ds;
    double A = init_param;
    double K = R*R+1.;
    double K_sq = K*K;
    double Kprime = 2.*K;
    return((A*cos(log(1. + (A*R)/exp((ds_sq)*(R_sq)))/C1)*(1. - 2.*(ds_sq)*(R_sq)))/
  ((exp((ds_sq)*(R_sq)) + A*R)*K) - 
 (C1*Kprime*sin(log(1. + (A*R)/exp((ds_sq)*(R_sq)))/C1))/(K_sq)); 
}

double initial_gauss_PHI2_m3bis_gen(double R,double init_param, vector<double> parameter) 
{
    double R_sq = R*R;
    double C2= parameter[0];
    double ds = parameter[1];
    double mu = parameter[2];
    double ds_sq = ds*ds;
    double A = init_param;
    double K = R*R+1.;
    double K_sq = K*K;
    double Kprime = 2.*K;
    return(-((C2*Kprime*cos(log(1. + (A*R)/exp((ds_sq)*(R_sq)))/C2))/(K_sq)) + 
 (A*(-1. + 2.*(ds_sq)*(R_sq))*sin(log(1. + (A*R)/exp((ds_sq)*(R_sq)))/C2))/
  ((exp((ds_sq)*(R_sq)) + A*R)*K)); 
}

double initial_parab_function(double x,double init_param, vector<double> parameter) 
{
    double a = -2*pow(x,2);
    return(a);
}

double line(double x,double t)
{
    return(4*t);
}


// initialize the fields at time zero
vector<vector<double>> initialize_fields(double d_min,double d_max,double dx,std::vector<double(*)(double,double,vector<double>)> funcs,double param_ic,vector<double> param,int gl, int gr,int ord)
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
            new_vect[n][j] = funcs[n](val,param_ic,param);
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
    fields_vect_new[j].insert(fields_vect_new[j].begin(),((-4*fields_vect_old[0][gl]+7.*fields_vect_old[0][gl+1]-4.*fields_vect_old[0][gl+2]+fields_vect_old[0][gl+3])/2/dx));
    //boundary conditions
    fields_vect_new[j].push_back((-(4*fields_vect_old[0][last_ind]-7.*fields_vect_old[0][last_ind-1]+4.*fields_vect_old[0][last_ind-2]-fields_vect_old[0][last_ind-3]) ));
    //fields_vect_new[j].push_back(((5.*fields_vect_old[j][last_ind]-11.*fields_vect_old[j][last_ind-1]+10.*fields_vect_old[j][last_ind-2]-5*    fields_vect_old[j][last_ind-3]+fields_vect_old[j][last_ind-4])/2/dx ));
}

void adv_boundaries_right(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //same of the bulk
    fields_vect_new[j][gl] = Dx[0](fields_vect_old[0],gl,dx) + artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt);
    //boundary conditions
    fields_vect_new[j][last_ind] = (-param[0]*Dx[0](fields_vect_old[0],last_ind,dx));
}

void adv_boundaries_left(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // boundary conditions
    fields_vect_new[j][gl] = (Dx[0](fields_vect_old[0],gl,dx)+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    //fields_vect_new[j].insert(fields_vect_new[j].begin(),(-3./2.*fields_vect_old[0][gl]+2.*fields_vect_old[0][gl+1]-1*2.*fields_vect_old[0][gl+2])+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) );
    
    

    // same of the bulk
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx));
    //fields_vect_new[j].push_back((-(3./2.*fields_vect_old[0][last_ind]-2.*fields_vect_old[0][last_ind-1]+1/2.*fields_vect_old[0][last_ind-2] )));
}



void adv_boundaries_periodic(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    //left boundary
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl+1]-fields_vect_old[0][last_ind])/2/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind-1])/2/dx ));
    
}

void adv_boundaries_periodic_backward_der(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[0].size()-1-gl;
    //left boundary
    
    fields_vect_new[j].insert(fields_vect_new[j].begin(),(-(fields_vect_old[0][gl]-fields_vect_old[0][last_ind])/dx+ artificial_diss(epsilon,ord,fields_vect_old,j,gl,dx,dt) ));   
    
    //right boundary
    fields_vect_new[j].push_back((-(fields_vect_old[0][last_ind]-fields_vect_old[0][last_ind-1])/dx ));

}


// wave equation

void refl_abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( (Dx[0](fields_vect_old[0],gl,dx)-Dx[0](fields_vect_old[1],gl,dx))/2 ));
    fields_vect_new[j].push_back((-(Dx[0](fields_vect_old[0],last_ind,dx)-Dx[0](fields_vect_old[1],last_ind,dx))/2));
}



void radiative_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
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

void radiative_outer_boundaries_PI_we(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gl;
    // at the left boundary, the origin, we don't put any conditions
    int ind_space_left = gl;
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space_left] = left_value;

    // right boundary conditions
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[j],last_ind,dx));
}


void abs_boundaries_PI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j].insert(fields_vect_new[j].begin(),( Dx[0](fields_vect_old[0],gl,dx)));
    fields_vect_new[j].push_back((-Dx[0](fields_vect_old[0],last_ind,dx)));
}

void abs_boundaries_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    fields_vect_new[j].insert(fields_vect_new[j].begin(), (Dx[0](fields_vect_old[1],gl,dx)));
    // right boundary conditions
    fields_vect_new[j].push_back(( -Dx[0](fields_vect_old[1],last_ind,dx)));
}

   
void no_boundary_conditions_PHI(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = left_value;
    //right no boundary
    fields_vect_new[j][last_ind] = right_value;   
}

void no_boundary_conditions_phi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = left_value;
    //right no boundary
    fields_vect_new[j][last_ind] = right_value;   
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
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);    
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
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);    
}

// R rescaling
void no_boundary_conditions_phi_hyp(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param)
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


void no_boundary_conditions_PI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}


void no_boundary_conditions_PHI_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
   
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_phi_hyp_Chi(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
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


void radiative_outer_boundaries_PI_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    
    
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    
    //double derivative = 3.* ( pow((x+dx),2)*fields_vect_old[1][ind_space+1] - pow((x-dx),2)*fields_vect_old[1][ind_space-1] )/(pow((x+dx),3)-pow((x-dx),3)) +(pow(fields_vect_old[1][ind_space],2)-pow(fields_vect_old[0][ind_space],2) );
    double left_value = evo(j,ind_space,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = ( left_value);

    
    
    // right boundary conditions
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx));

    
}

//------------------- MODEL 1: CHARACTERISTIC VARIABLE, FLAT SLICES -------------------//
//fields_vect [0]-> PsiPlus, [1]->PsiMinus, [2]->Psi

void radiative_outer_boundaries_PsiPlus_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = 0;  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void radiative_outer_boundaries_PsiMinus_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = -Dx[0](fields_vect_old[0],last_ind,dx) - Dx[0](fields_vect_old[1],last_ind,dx);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}
//------------------- MODEL 1: CHARACTERISTIC VARIABLE, CHI RESCALING, HYPERBOLOIDAL FOLIATION -------------------//
//fields_vect [0]-> PsiPlus, [1]->PsiMinus, [2]->Psi

void no_boundary_conditions_PsiPlus_charvar_hyp_Chi_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_PsiMinus_charvar_hyp_Chi_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi_charvar_hyp_Chi_m1(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

// ------------- // MODEL 2 // ------------- //


void radiative_outer_boundaries_PI1_m2(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx) );
}

void radiative_outer_boundaries_PI2_m2(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
                        
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[3],last_ind,dx) );
}

// ------------- // MODEL 3 // ------------- //


void radiative_outer_boundaries_PI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)-fields_vect_old[0][last_ind]/x );
}

void radiative_outer_boundaries_PI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
                        
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[3],last_ind,dx)-fields_vect_old[3][last_ind]/x );
}


void no_boundary_conditions_PHI1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{    
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
   
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_PHI2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
   
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}


void no_boundary_conditions_phi1_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
   
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_phi2_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
   
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  

    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

// ------------- // MODEL 3bis Generalized // ------------- //

void radiative_outer_boundaries_PI1_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[0],last_ind,dx)-fields_vect_old[0][last_ind]/x );
}

void radiative_outer_boundaries_PI2_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    // at the left boundary, the origin, we don't put any conditions
    int ind_space = gl;
    double x = dmin;
    double K = x*x+1.;
    double dK = 2.*x;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;

    
    // right boundary conditions
    x = dmax;
   
    int last_ind = fields_vect_old[j].size()-1-gr;
    fields_vect_new[j][last_ind] = (-Dx[0](fields_vect_old[3],last_ind,dx)-fields_vect_old[3][last_ind]/x );
}


void no_boundary_conditions_PHI1_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;
    // right boundary conditions
    x = dmax;
    
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl); 
    fields_vect_new[j][last_ind] = right_value;
}

void no_boundary_conditions_PHI2_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;  
    // left boundary conditions
    int ind_space = gl;
    double x = dmin;

    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;
    
    // right boundary conditions
    x = dmax;
   
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl); 
    fields_vect_new[j][last_ind] = right_value;
}


void no_boundary_conditions_phi1_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    int ind_space = gl;
    double x = dmin;
    
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;
    // right boundary conditions
    x = dmax;

    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl); 
    fields_vect_new[j][last_ind] = right_value;
}

void no_boundary_conditions_phi2_m3bis_gen(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,std::vector<double (*)(std::vector<double>,int,double)> Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    // left boundary conditions
    int ind_space = gl;
    double x = dmin;
    double derivative = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    fields_vect_new[j][ind_space] = derivative;
    // right boundary conditions
    x = dmax;

    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl); 
    fields_vect_new[j][last_ind] = right_value;
}

// -----------  HYPERBOLOIDAL COMPACTIFIED CHI RESCALED M3 NO BOUNDARY CONDITIONS ----------- //

void no_boundary_conditions_Psi1_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}
void no_boundary_conditions_Psi2_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi1Plus_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi1Minus_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi2Plus_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

void no_boundary_conditions_Psi2Minus_charvar_hyp_Chi_m3(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double s = param[0]; 
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
}

// Curved wave equation //
void no_boundary_conditions_F(std::vector<std::vector<double>> &fields_vect_new,std::vector<std::vector<double>> fields_vect_old,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,derivative_vector Dx,artificial_dissipation_function artificial_diss,double epsilon,int ord,std::vector<double> &param, evolution_function evo)
{
    int last_ind = fields_vect_old[j].size()-1-gr;
    int ind_space_left=gl;
    double left_value = evo(j,gl,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);
    double right_value = evo(j,last_ind,fields_vect_old,dx,dmin,param,t,Dx,artificial_diss,epsilon,ord,dt,gl);  
    //left no boundary
    fields_vect_new[j][ind_space_left] = (left_value);
    //right no boundary
    fields_vect_new[j][last_ind] = (right_value);   
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

void ghost_point_extrapolation_1_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation1(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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


void ghost_point_extrapolation2_TEM(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_2_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int last_ind = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    for (int i=gl-1; i>=0;i--)
    {
        field_vect[j][i] = field_vect[j][i+1]-(-3./2.*field_vect[j][i+1]+2.*field_vect[j][i+2]-1./2.*field_vect[j][i+3])+1./2.*(2.*field_vect[j][i+1]-5.*field_vect[j][i+2]+4.*field_vect[j][i+3]-field_vect[j][i+4]);
    }
        
    for (int i=last_ind-gr; i<last_ind;i++)
    {
        //cout<<i<<endl;
        field_vect[j][i] = field_vect[j][i-1]+(3./2.*field_vect[j][i-1]-2.*field_vect[j][i-2]+1./2.*field_vect[j][i-3])+1./2.*(2.*field_vect[j][i-1]-5.*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4]);
    }
}

void ghost_point_extrapolation_4_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_6_ord(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation3(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_1_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_4_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2 || j==3 || j==5)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1 || j==4)
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

void ghost_point_extrapolation_1_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2 || j==3 || j==5)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1 || j==4)
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

void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry_m3bis(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==1 || j==0 ||  j==3|| j==5 )
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==2 ||j==4)
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
void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry_curved_element(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int S = field_vect[j].size();
    // attachment of the ghost points to the boundary of the function
    // PI is even
    if (j==0 || j==2 ||j==3)
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
    /*
    if (j==3)
    {
        for (int i=gl-1; i>=0;i--)
        {
            field_vect[j][i] = field_vect[j][i+1]-(-3./2.*field_vect[j][i+1]+2.*field_vect[j][i+2]-1./2.*field_vect[j][i+3])+1./2.*(2.*field_vect[j][i+1]-5.*field_vect[j][i+2]+4.*field_vect[j][i+3]-field_vect[j][i+4]);
        }
    }*/
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = 4.*field_vect[j][i-1]-6*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4];
    }
}


void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry_charvar(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    
    // attachment of the ghost points to the boundary of the function
    // PsI is even
    if (j==2 )
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    //PsiPlus condition
    if (j==0 )
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[0][i] = field_vect[1][2*gl-i];
        }
        
    }
    // PsiMinus condition
    if (j==1 )
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[1][i] = field_vect[0][2*gl-i];
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = 4.*field_vect[j][i-1]-6*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4];
    }
}

void ghost_point_extrapolation_2_ord_TEM_spherical_symmetry_charvar_Chi(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
{
    int S = field_vect[j].size();
    double s = param[0];
    double r;
    double Chi;
    double R;
    // attachment of the ghost points to the boundary of the function
    // PsI is even
    if (j==2 || j==5 )
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    //PsiPlus condition
    if (j==0 || j==3 )
    {
        for (int i=0; i<gl;i++)
        {
            r = dmin+dx*double(gl-i);
            R = r/(1.-r*r/s/s);
            Chi = sqrt(1.+R*R);
            field_vect[j][i] = Chi*field_vect[j+1][2*gl-i];
        }
    }
    // PsiMinus condition
    if (j==1 || j==4)
    {
        for (int i=0; i<gl;i++)
        {
            r = dmin+dx*double(gl-i);
            R = r/(1.-r*r/s/s);
            Chi = sqrt(1.+R*R);
            field_vect[j][i] = field_vect[j-1][2*gl-i]/Chi;
        }
    }
    
    for (int i=field_vect[j].size()-gr; i<field_vect[j].size();i++)
    {
        field_vect[j][i] = 4.*field_vect[j][i-1]-6*field_vect[j][i-2]+4.*field_vect[j][i-3]-field_vect[j][i-4];
    }
    
}

void ghost_point_extrapolation_2_ord_spline_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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
    if (j==0 || j==2 || j==3 || j==5)
    {
        for (int i=0; i<gl;i++)
        {
            field_vect[j][i] = field_vect[j][2*gl-i];
        }
    }
    
    // PHI is odd
    if (j==1 || j==4)
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

void ghost_point_extrapolation_4_ord_spherical_symmetry_rescaled(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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

void ghost_point_extrapolation_6_ord_spherical_symmetry(std::vector<std::vector<double>> &field_vect,double t,double dx, double dt, int j,int gl, int gr,double dmin,double dmax,vector<double> param)
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
double artificial_dissipation_2(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*(dt/dx)*pow(-1,ord)*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_Husa(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon/pow(dx,1)*pow(-1,ord)/16*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
}

// to control!
double artificial_dissipation_4(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*pow(dx,3)*pow(-1,ord)/4*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4*copy_fields_vect[j][i-1]) );
}

double artificial_dissipation_2_ani(double epsilon,int ord,std::vector<std::vector<double>> copy_fields_vect,int j,int i,double dx,double dt)      
{
    return(-epsilon*pow(dx,-1.)*pow(-1,ord)/4.*(copy_fields_vect[j][i+2]-4.*copy_fields_vect[j][i+1]+6.*copy_fields_vect[j][i]+copy_fields_vect[j][i-2]-4.*copy_fields_vect[j][i-1]) );
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



void print_f(std::vector< std::vector<double> > fields_vect, double dmin, double dmax, double dx, string name_file,string name_folder, int gl, int gr,MPI_Status status, int totalnodes, int mynode,MPI_Request request)
{
    int S =  int( ((dmax+dx*double(gr))-(dmin-dx*double(gl)) +dx/2 ) / dx) + 1;
    double index_dmin_local,index_dmax_local;
        
    // defining the first and last spatial index of the domain, the actual process will work inside this range
        
    index_dmin_local = (gl+1) + mynode * int((S-2-gl-gr)/totalnodes);
    if (mynode==totalnodes-1)
    {
        index_dmax_local =  S-1-gr;
    }
    else
    {
        index_dmax_local = (gl+1) + (mynode+1) * int((S-2-gl-gr)/totalnodes);
    }
    ofstream myfile_print_f;
    myfile_print_f.open (name_file,ios::app);
    myfile_print_f.precision(20);
    // headers of the columns
    myfile_print_f<<"x,";
    for (int i=0; i<fields_vect.size()-1; i++)
    {
        myfile_print_f << "field"<<to_string(i)<<",";
    }    
    myfile_print_f << "field"<<to_string(fields_vect.size()-1)<<"\n";
    
    // for every spatial point
    if(mynode==0)
    {
        for(int j=0;j<gl+1;j++)
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


// Functions usefull in a parallelized setting

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


// ---------- Model setting function ---------- //

// INITIAL CONDITIONS //

// ---  wave equation --- //


void initial_conditions_setting_we_null_field(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI
    initial_conditions_vector.push_back(&initial_null);        //PHI       
    initial_conditions_vector.push_back(&initial_null);        //phi
}



void initial_conditions_setting_we_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_null);    //PI
    initial_conditions_vector.push_back(&initial_gauss_PHI);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss);        //phi
}
    
// ---  model 1 --- //

    
void initial_conditions_setting_m1(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss_m1);    //PI
    initial_conditions_vector.push_back(&initial_null);        //PHI       
    initial_conditions_vector.push_back(&initial_null);        //phi
}

// ---  model 2 --- //

void initial_conditions_setting_m2_null_field(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);//PI1
    initial_conditions_vector.push_back(&initial_null);        //PHI1      
    initial_conditions_vector.push_back(&initial_null);        //phi1
    initial_conditions_vector.push_back(&initial_gauss_m2);    //PI2
    initial_conditions_vector.push_back(&initial_null);        //PHI2       
    initial_conditions_vector.push_back(&initial_null);        //phi2
}


void initial_conditions_setting_m2_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_null);//PI1
    initial_conditions_vector.push_back(&initial_gauss_PHI);        //PHI1      
    initial_conditions_vector.push_back(&initial_gauss);        //phi1
    initial_conditions_vector.push_back(&initial_null);    //PI2
    initial_conditions_vector.push_back(&initial_gauss_m2_PHI);        //PHI2       
    initial_conditions_vector.push_back(&initial_gauss_m2);        //phi2
}

// ---  model 3 --- //


void initial_conditions_setting_m3_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_null);//PI1
    initial_conditions_vector.push_back(&initial_gauss_PHI);        //PHI1      
    initial_conditions_vector.push_back(&initial_gauss);        //phi1
    initial_conditions_vector.push_back(&initial_null);    //PI2
    initial_conditions_vector.push_back(&initial_gauss_m3_PHI);        //PHI2       
    initial_conditions_vector.push_back(&initial_gauss_m3);        //phi2
}


void initial_conditions_setting_m3(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI1
    initial_conditions_vector.push_back(&initial_null);        //PHI1       
    initial_conditions_vector.push_back(&initial_null);        //phi1
    initial_conditions_vector.push_back(&initial_null);    //PI2
    initial_conditions_vector.push_back(&initial_null);        //PHI2
    initial_conditions_vector.push_back(&initial_unity);        //phi2
}

// ---  model 4 --- //


void initial_conditions_setting_m4_null_field(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI1
    initial_conditions_vector.push_back(&initial_null);        //PHI1       
    initial_conditions_vector.push_back(&initial_null);        //phi1
    initial_conditions_vector.push_back(&initial_null);    //PI2
    initial_conditions_vector.push_back(&initial_null);        //PHI2
    initial_conditions_vector.push_back(&initial_null);        //phi2
}

// --- hyperboloidal wave equation --- //

void initial_conditions_setting_hyp_Chi_we_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss_PI_compactified_Chi_we_null_time_derivative);    //PI
    initial_conditions_vector.push_back(&initial_gauss_PHI_compactified_Chi_we_null_time_derivative);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss_phi_compactified_Chi_we_null_time_derivative);        //phi
}

// --- hyperboloidal model 1 --- //

void initial_conditions_setting_hyp_Chi_m1_charvar_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss_PsiPlus_compactified_Chi_charvar_m1_null_field);    //PI
    initial_conditions_vector.push_back(&initial_gauss_PsiMinus_compactified_Chi_charvar_m1_null_field);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss_Psi_compactified_Chi_charvar_m1_null_field);        //phi
}

// --- hyperboloidal model 3 --- //

void initial_conditions_setting_hyp_Chi_m3_charvar_null_time_derivative(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss_Psi1Plus_compactified_Chi_charvar_m3_null_field);    //PI
    initial_conditions_vector.push_back(&initial_gauss_Psi1Minus_compactified_Chi_charvar_m3_null_field);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss_Psi1_compactified_Chi_charvar_m3_null_field); 
    initial_conditions_vector.push_back(&initial_gauss_Psi2Plus_compactified_Chi_charvar_m3_null_field);    //PI
    initial_conditions_vector.push_back(&initial_gauss_Psi2Minus_compactified_Chi_charvar_m3_null_field);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss_Psi2_compactified_Chi_charvar_m3_null_field); //phi
}


// mofidified models
void initial_conditions_setting_hyp_Chi_we_null_field(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss_PI_compactified_Chi_we_null_field);    //PI
    initial_conditions_vector.push_back(&initial_gauss_PHI_compactified_Chi_we_null_field);        //PHI       
    initial_conditions_vector.push_back(&initial_gauss_phi_compactified_Chi_we_null_field);        //phi
}

void initial_conditions_setting_curved_we_simple(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI
    initial_conditions_vector.push_back(&initial_null);        //PHI    
    initial_conditions_vector.push_back(&initial_null);        //phi
    initial_conditions_vector.push_back(&initial_unity);
}



void initial_conditions_setting_m3hyb(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI1
    initial_conditions_vector.push_back(&initial_null);        //PHI1       
    initial_conditions_vector.push_back(&initial_null);        //phi1
    initial_conditions_vector.push_back(&initial_gauss);    //PI2
    initial_conditions_vector.push_back(&initial_null);        //PHI2
    initial_conditions_vector.push_back(&initial_null);        //phi2
}

void initial_conditions_setting_DGModel(std::vector<double(*)(double, double,vector<double>)> &initial_conditions_vector)
{
    initial_conditions_vector.push_back(&initial_gauss);    //PI1
    initial_conditions_vector.push_back(&initial_null);        //PHI1       
    initial_conditions_vector.push_back(&initial_null);        //phi1
    initial_conditions_vector.push_back(&initial_gauss);    //PI2
    initial_conditions_vector.push_back(&initial_null);        //PHI2
    initial_conditions_vector.push_back(&initial_null);        //phi2
}

// EVOLUTION VECTOR //
void evo_vector_setting_we(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&wave_eq_spherical_PI);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    
}

void evo_vector_setting_m1(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model1_PI);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi); 
}

void evo_vector_setting_m2(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model2_PI1);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    evo_vector.push_back(&model2_PI2);
    evo_vector.push_back(&model3_PHI2);
    evo_vector.push_back(&model3_phi2);
}

void evo_vector_setting_m1bis(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model1_PI_bis);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    
}

void evo_vector_setting_m3(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model3_PI1);
    evo_vector.push_back(&model3_PHI1);
    evo_vector.push_back(&model3_phi1);
    evo_vector.push_back(&model3_PI2);
    evo_vector.push_back(&model3_PHI2);
    evo_vector.push_back(&model3_phi2);
}



void evo_vector_setting_hyp_we_Chi(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&wave_eq_compactified_PI_Chi);
    evo_vector.push_back(&wave_eq_compactified_PHI_Chi);
    evo_vector.push_back(&wave_eq_compactified_phi_Chi);
    
}

void evo_vector_setting_hyp_m1_Chi_charvar(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model1_charvar_compactified_PsiPlus_Chi);
    evo_vector.push_back(&model1_charvar_compactified_PsiMinus_Chi);
    evo_vector.push_back(&model1_charvar_compactified_Psi_Chi);
    
}

void evo_vector_setting_hyp_m3_Chi_charvar(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model3_charvar_compactified_PsiOnePlus_Chi);
    evo_vector.push_back(&model3_charvar_compactified_PsiOneMinus_Chi);
    evo_vector.push_back(&model3_charvar_compactified_PsiOne_Chi);
    evo_vector.push_back(&model3_charvar_compactified_PsiTwoPlus_Chi);
    evo_vector.push_back(&model3_charvar_compactified_PsiTwoMinus_Chi);
    evo_vector.push_back(&model3_charvar_compactified_PsiTwo_Chi);
    
}



void evo_vector_setting_curved_we_simple(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&curved_wave_eq_spherical_PI);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    evo_vector.push_back(&curved_wave_eq_spherical_F_simple);
    
}



void evo_vector_setting_curved_we_simple_not_singular(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&curved_wave_eq_spherical_PI_1);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    evo_vector.push_back(&curved_wave_eq_spherical_F_simple);
    
}

void evo_vector_setting_DGModel(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&wave_eq_spherical_PI);
    evo_vector.push_back(&wave_eq_PHI);
    evo_vector.push_back(&wave_eq_phi);
    evo_vector.push_back(&DG_Model2_PI2);
    evo_vector.push_back(&model3_PHI2);
    evo_vector.push_back(&model3_phi2);
    
}


void evo_vector_setting_m3Hyb(std::vector< evolution_function > &evo_vector)
{
    evo_vector.push_back(&model3Hyb_PI1);
    evo_vector.push_back(&model3_PHI1);
    evo_vector.push_back(&model3_phi1);
    evo_vector.push_back(&model3Hyb_PI2);
    evo_vector.push_back(&model3_PHI2);
    evo_vector.push_back(&model3_phi2);
}

// BOUNDARY CONDITIONS //
void bc_vector_setting_we(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI1_m3);
    b_func_vec[1] = (no_boundary_conditions_PHI);
    b_func_vec[2] = (no_boundary_conditions_phi);
}

void bc_vector_setting_m1(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI_m1);
    b_func_vec[1] = (no_boundary_conditions_PHI);
    b_func_vec[2] = (no_boundary_conditions_phi);
}

void bc_vector_setting_m2(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI1_m2);
    b_func_vec[1] = (no_boundary_conditions_PHI1_m3);
    b_func_vec[2] = (no_boundary_conditions_phi1_m3);
    b_func_vec[3] = (radiative_outer_boundaries_PI2_m2);
    b_func_vec[4] = (no_boundary_conditions_PHI2_m3);
    b_func_vec[5] = (no_boundary_conditions_phi2_m3);
}

void bc_vector_setting_m3(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI1_m3);
    b_func_vec[1] = (no_boundary_conditions_PHI1_m3);
    b_func_vec[2] = (no_boundary_conditions_phi1_m3);
    b_func_vec[3] = (radiative_outer_boundaries_PI2_m3);
    b_func_vec[4] = (no_boundary_conditions_PHI2_m3);
    b_func_vec[5] = (no_boundary_conditions_phi2_m3);
}

void bc_vector_setting_hyp_we_Chi(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (no_boundary_conditions_PI_hyp_Chi);
    b_func_vec[1] = (no_boundary_conditions_PHI_hyp_Chi);
    b_func_vec[2] = (no_boundary_conditions_phi_hyp_Chi);
}

void bc_vector_setting_hyp_m1_Chi_charvar(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (no_boundary_conditions_PsiPlus_charvar_hyp_Chi_m1);
    b_func_vec[1] = (no_boundary_conditions_PsiMinus_charvar_hyp_Chi_m1);
    b_func_vec[2] = (no_boundary_conditions_Psi_charvar_hyp_Chi_m1);
}

void bc_vector_setting_hyp_m3_Chi_charvar(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (no_boundary_conditions_Psi1Plus_charvar_hyp_Chi_m3);
    b_func_vec[1] = (no_boundary_conditions_Psi1Minus_charvar_hyp_Chi_m3);
    b_func_vec[2] = (no_boundary_conditions_Psi1_charvar_hyp_Chi_m3);
    b_func_vec[3] = (no_boundary_conditions_Psi2Plus_charvar_hyp_Chi_m3);
    b_func_vec[4] = (no_boundary_conditions_Psi2Minus_charvar_hyp_Chi_m3);
    b_func_vec[5] = (no_boundary_conditions_Psi2_charvar_hyp_Chi_m3);
}



// ----------- modified models ------------//

void bc_vector_setting_hyp_curved_we_simple(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI_we);
    b_func_vec[1] = (no_boundary_conditions_PHI);
    b_func_vec[2] = (no_boundary_conditions_phi);
    b_func_vec[3] = (no_boundary_conditions_F);
}


void bc_vector_setting_DGModel(std::vector< boundary_conditions_function> &b_func_vec)
{
    b_func_vec[0] = (radiative_outer_boundaries_PI1_m3);
    b_func_vec[1] = (no_boundary_conditions_PHI1_m3);
    b_func_vec[2] = (no_boundary_conditions_phi1_m3);
    b_func_vec[3] = (radiative_outer_boundaries_PI2_m3);
    b_func_vec[4] = (no_boundary_conditions_PHI2_m3);
    b_func_vec[5] = (no_boundary_conditions_phi2_m3);
    
}
