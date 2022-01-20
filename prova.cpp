// Epa! //
#include <iostream>
#include <fstream>
#include<vector>
#include <math.h>
#include <chrono>
#include <mpi.h>
#include <chrono>

#include "boost/multi_array.hpp"
using namespace std;
using namespace std::chrono;

typedef vector<vector<vector< vector<double> >>> grid;
typedef vector<vector<vector< vector<double> >>> fields_vector;
typedef double(*partial_derivative_operator)(fields_vector fields_vec,int ind_field,int i, int j, int k,double dx);
typedef std::vector<partial_derivative_operator>  derivatives_vector;
typedef double(*evolution_function)(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );
typedef void (*communication_function)(fields_vector &fields_vect,grid Grid,int n,vector<int> index_dmax_local,vector<int> index_dmin_local,int mynode,int totalnodes,MPI_Status status);
typedef void(*ghost_point_extrapolation_function)(fields_vector &field_vect,grid Grid, vector<int> S, double t,vector<double> dx, double dt,vector<int> gl, vector<int> gr,vector<double> dmin,vector<double> dmax,vector<int> index_dmin_local, vector<int> index_dmax_local,int mynode, int totalnodes);
const int dim = 3;

double partial_dx(fields_vector fields_vec,int ind_field,int i, int j, int k,double dx);

double partial_dy(fields_vector fields_vec,int ind_field,int i, int j, int k,double dy);

double partial_dz(fields_vector fields_vec,int ind_field,int i, int j, int k,double dz);


grid initialization_grid(vector<double> dmin,vector<double> dmax,vector<double> dx,vector<int> gl,vector<int> gr,vector<int> ord,int dim);

grid initialize_fields(vector<double> dmin,vector<double> dmax,vector<double> dx,vector<int> gl, vector<int> gr,vector<int> ord,int dim,grid Grid,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> param_ic);

double we_PI_2D(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double we_PI_3D(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double we_PHIx(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double we_PHIy(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double we_PHIz(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double we_phi(int field_ind,int i, int j, int k,fields_vector fields_vect, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double null_eq(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt );

double initial_gauss(vector<double> x,vector<double> init_param);

double initial_gauss_PHIx(vector<double> x,vector<double> init_param);

double initial_gauss_PHIy(vector<double> x,vector<double> init_param);

double initial_null(vector<double> x,vector<double> init_param);

void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid, vector<int> S, double t,vector<double> dx, double dt,vector<int> gl, vector<int> gr,vector<double> dmin,vector<double> dmax,vector<int> index_dmin_local, vector<int> index_dmax_local,int mynode, int totalnodes);

int print_f(fields_vector fields_vect, grid Grid, vector<double> dmin, vector<double> dmax, vector<double> dx, string name_file,string name_folder, vector<int> gl, vector<int> gr,int mynode, int totalnodes);

fields_vector onestep_RK4_1(fields_vector fields_vect,grid Grid,vector<double> dmin,vector<double> dmax,vector<double> dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,double t, vector<int> gl, vector<int> gr, ghost_point_extrapolation_function ghost_point_extrapolation,derivatives_vector Dx,vector<int> ord,communication_function communication, MPI_Status status, int totalnodes, int mynode);
//,std::vector< boundary_conditions_function > bc, ghost_point_extrapolation_function ghost_point_extrapolation, artificial_dissipation_function artificial_diss,double epsilon, communication_function communication

void communication(fields_vector &fields_vect,grid Grid,int n,vector<int> index_dmax_local,vector<int> index_dmin_local,int mynode,int totalnodes,MPI_Status status);


void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid,double t,vector<double> dx, double dt, int j,vector<int> gl, vector<int> gr,vector<double> dmin,vector<double> dmax);

// ------------------main function ------------------ //

int main(int argc, char **argv) 
{
     auto start = high_resolution_clock::now();
    // set the decimal numbers precision of the output
     cout.precision(10);
    int mynode, totalnodes;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    
   
    
    
    vector<double> dmin {-1.5,-1.5,0};
    vector<double> dmax {1.5,1.5,0};
    vector<double> dx {0.1,0.1,0};
    vector<int> gl = {2,2,0}, gr = {2,2,0}, ord = {2,2,0};
    int step_to_save = 30;
    double dt = 0.05;
    double integration_interval = dt*step_to_save;
    
    vector<double> ic_par = {1};
    vector<double> param = {1};
    
    string file_path = "./data/prova/data0/";
    string name_file = file_path + "processor_"+to_string(mynode)+"_ampl_"+to_string(ic_par[0])+"_eps"+to_string(0)+"_dx_"+to_string(dx[0])+"steps"+to_string(step_to_save)+"last_time"+to_string(integration_interval)+".csv";
    
    vector<double(*)(vector<double>, vector<double>)> initial_conditions;
    initial_conditions.push_back(initial_gauss_PHIx);    
    initial_conditions.push_back(initial_gauss_PHIy);        
    initial_conditions.push_back(initial_null);       
    initial_conditions.push_back(initial_null); 
    initial_conditions.push_back(initial_gauss);
    
    vector<evolution_function> R_vect;
    R_vect.push_back(we_PHIx);
    R_vect.push_back(we_PHIy);
    R_vect.push_back(null_eq);
    R_vect.push_back(we_PI_2D);
    R_vect.push_back(we_phi);
    
    derivatives_vector Dx;
    Dx.push_back(partial_dx);
    Dx.push_back(partial_dy);
    Dx.push_back(partial_dz);
    grid Grid= initialization_grid(dmin,dmax,dx, gl,gr,ord,dim);
    //cout<<"x_0 size: "<<Grid.size()<<endl;
    //cout<<"x_1 size: "<<Grid[0].size()<<endl;
    //cout<<"x_2 size: "<<Grid[0][0].size()<<endl;
    //cout<<"dim size: "<<Grid[0][0][0].size()<<endl;
    fields_vector new_fields_vector = initialize_fields(dmin,dmax,dx, gl,  gr, ord, dim, Grid,initial_conditions,ic_par);
    
    int N_points;
    
    for(double t =0;t<integration_interval;t=t+dt)
    {
        cout<<"time :"<<t<<endl;
        N_points = print_f(new_fields_vector, Grid, dmin, dmax, dx, name_file, name_file,  gl,  gr,mynode,totalnodes);
        cout<<"N point "<<(N_points)<<endl;
        new_fields_vector = onestep_RK4_1(new_fields_vector, Grid,dmin,dmax,dx,param,dt,R_vect,t,gl,gr,ghost_point_extrapolation_4_ord,Dx,ord,communication,status,totalnodes,mynode);
        
    }
    ofstream myfile_print_f;
    myfile_print_f.open (name_file,ios::app);
    myfile_print_f<<N_points;
    myfile_print_f.close();
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by processor "<<mynode<<" : "<< duration.count()/1000000. << " seconds" << endl;
    MPI_Finalize();
    return 0;
} 



// ----------------------- FUNCIONS DEFINITIONS ----------------------- //





fields_vector onestep_RK4_1(fields_vector fields_vect,grid Grid,vector<double> dmin,vector<double> dmax,vector<double> dx,std::vector<double> param, double dt, std::vector< evolution_function > evo,double t, vector<int> gl, vector<int> gr,ghost_point_extrapolation_function ghost_point_extrapolation,derivatives_vector Dx,vector<int> ord,communication_function communication,MPI_Status status, int totalnodes, int mynode)  
{
    // setting the ghost point equal to the order of the artificial dissipator
    for(int d=0;d<dim;d++)
    {
        if (ord[d]>gl[d])
        {
            gl[d] = ord[d];
        }
        if(ord[d]>gr[d])
        {
            gr[d] = ord[d];
        }
    }
    //cout<<" processor: "<<mynode<<"time "<<t<<endl;
    int nitems=2;
    int N = fields_vect.size();
    int dim_x = Grid.size(),dim_y=Grid[0].size(),dim_z=Grid[0][0].size();

    
    fields_vector k1 (N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector support_k1(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector k2(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector support_k2(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector k3(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector support_k3(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector k4(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
    fields_vector support_k4(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));
       
        
    // defining the first and last spatial index of the domain, the actual process will work inside this range
    // inside the core domain there are also the borders, but the GP are left out of this domain
        
    vector<int> S(dim);
    vector<int> index_dmin_local(dim),index_dmax_local(dim);
    for(int d=0;d<dim;d++)
    {
        if(dmax[d]!=dmin[d]) // we must consider the case of a dimension be suppressed
        {
            S[d] =  int( ((dmax[d]+dx[d]*double(gr[d]))-(dmin[d]-dx[d]*double(gl[d])) +dx[d]/2 ) / dx[d]) + 1;
        }
        else
        {
            S[d] = 1;
        }
        //cout<<"S"<<d<<": "<<S[d]<<endl;
        // defining the first and last spatial index of the domain, the actual process will work inside this range
        if (S[d]!=1)
        {
            index_dmin_local[d] = (gl[d]) + mynode * int((S[d]-gl[d]-gr[d])/totalnodes);
            //cout<<"mx from processor: "<<mynode<<" dmin local: dim "<<d<<" index "<<index_dmin_local[d]<<endl;
            if (mynode==totalnodes-1)
            {
                index_dmax_local[d] =  S[d]-gr[d];
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
            }
            else
            {
                index_dmax_local[d] = (gl[d]) + (mynode+1) * int((S[d]-gl[d]-gr[d])/totalnodes);
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
            }
        }
        else
        {
            index_dmin_local[d] = 0;
            index_dmax_local[d] = 1;
        }
    }
    
    
    // ! we parallelize on the first coordinate !
    for(int d=1;d<dim;d++)
    {
        if (S[d]!=1)
        {
            index_dmin_local[d] = (gl[d]) ;
            //cout<<"mx from processor: "<<mynode<<" dmin local: dim "<<d<<" index "<<index_dmin_local[d]<<endl;
           
            index_dmax_local[d] =  S[d]-gr[d];
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
        }
        else
        {
            index_dmin_local[d] = 0;
            index_dmax_local[d] = 1;
        }
    }

    
    // populate the ghost zone of the fields, we need them to calculate k1
    /*
    for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }*/
    
    // k1 building
    for (int n=0; n <N; n++)
    {
        // evualuating the "bulk" of k1
        // we have to consider the physical domain of the fields (so exclude the GP) and then exclude the boundaries value (+1 and -1)
        // we are also dividing the spatial domain in subdomains, in which each processor run in parallel
        
        // Bulk processors evolving routine for k1 (not optimized, it's computed also on the borders but it is then overwritten with thr bc func
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    k1[n][i][j][k] = evo[n](n,i,j,k,fields_vect,Grid,dx,param,t,Dx,dt);
                    // !!! artificial dissipation !!!
                }
            }
        }
        
        // evaluating the boundary of k1
        //-------------bc[n](k1,fields_vect,Grid,t,dx,dt,n,gl,gr,dmin,dmax,Dx,param,evo[j]);
        //k1[j][gl] += artificial_diss(epsilon,ord,fields_vect,j,gl,dx,dt);
        //k1[j][gr] += artificial_diss(epsilon,ord,fields_vect,j,gr,dx,dt);
        
        // computing the argument for the next coefficient k2
        
       
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    support_k1[n][i][j][k] =  (k1[n][i][j][k])*dt/2. + fields_vect[n][i][j][k];
                }
            }
        }
        
    
    
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(support_k1,Grid,n,index_dmax_local,index_dmin_local,mynode,totalnodes,status);
        // --------ghost_point_extrapolation(k1, t,dx,dt,j,gl,gr,dmin,dmax);
    }
    
    // GP extrapolation for support_k1 vector
    ghost_point_extrapolation(support_k1,Grid,S,t,dx,dt,gl,gr,dmin,dmax,index_dmin_local,index_dmax_local,mynode,totalnodes);
    
    
    // k2 building
    for (int n=0; n <N; n++)
    {
        // Bulk processors evolving routine for k1
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)  // including the boundaries 
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    k2[n][i][j][k] = evo[n](n,i,j,k,support_k1,Grid,dx,param,t+dt/2,Dx,dt/2);
                    //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2));
                }
            }
        }
        
        
       
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    support_k2[n][i][j][k] = (k2[n][i][j][k])*dt/2. + fields_vect[n][i][j][k];
                }
            }
        }
        

        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(support_k2,Grid,n,index_dmax_local,index_dmin_local,mynode,totalnodes,status);
        //boundary conditions update
        // ---------bc[j](k2,support_k1, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        // artificial dissipation
        //k2[j][gl] += artificial_diss(epsilon,ord,support_k1,j,gl,dx,dt/2);
        //k2[j][gr] += artificial_diss(epsilon,ord,support_k1,j,gr,dx,dt/2);
        
        
        // ------------- ghost_point_extrapolation(k2, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax);
        
    }
    ghost_point_extrapolation(support_k2,Grid,S,t,dx,dt,gl,gr,dmin,dmax,index_dmin_local,index_dmax_local,mynode,totalnodes);
    
    // k3 building
    for (int n=0; n <N; n++)
    {
        // Bulk processors evolving routine for k1
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)  // including the boundaries 
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    k3[n][i][j][k] = evo[n](n,i,j,k,support_k2,Grid,dx,param,t+dt/2,Dx,dt/2);
                    //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2));
                }
            }
        }
        
        
        //----bc[j](k3,support_k2, t+dt/2.,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        //k3[j][gl] += artificial_diss(epsilon,ord,support_k2,j,gl,dx,dt/2);
        //k3[j][gr] += artificial_diss(epsilon,ord,support_k2,j,gr,dx,dt/2);
        
        //------ghost_point_extrapolation(k3, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax);
        
       
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    support_k3[n][i][j][k] = k3[n][i][j][k]*dt + fields_vect[n][i][j][k];
                }
            }
        }        
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        
        communication(support_k3,Grid,n,index_dmax_local,index_dmin_local,mynode,totalnodes,status);
        //------ghost_point_extrapolation(k3, t+dt/2.,dx,dt/2,j,gl,gr,dmin,dmax);
    }
   ghost_point_extrapolation(support_k3,Grid,S,t,dx,dt,gl,gr,dmin,dmax,index_dmin_local,index_dmax_local,mynode,totalnodes);
   
    // k4 building    
    for (int n=0; n <N; n++)
    {
        // Bulk processors evolving routine for k1
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)  // including the boundaries 
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    k4[n][i][j][k]=evo[n](n,i,j,k,support_k3,Grid,dx,param,t+dt,Dx,dt);
                    //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt/2));
                }
            }
        }
        
        // processes COMMUNICATION blocks // since we usually use a centered second order finite difference scheme, 
        // each processor needs to receive (and send) the borders of its subdomain
        communication(k4,Grid,n,index_dmax_local,index_dmin_local,mynode,totalnodes,status);
        //--------bc[j](k4,support_k3, t+dt,dx,dt,j,gl,gr,dmin,dmax,Dx,artificial_diss,epsilon,ord,param,evo[j]);
        //k4[j][gl] += artificial_diss(epsilon,ord,support_k3,j,gl,dx,dt);
        //k4[j][gr] += artificial_diss(epsilon,ord,support_k3,j,gr,dx,dt);
        
        //-----ghost_point_extrapolation(k4, t+dt,dx,dt,j,gl,gr,dmin,dmax);
    }
    ghost_point_extrapolation(k4,Grid,S,t,dx,dt,gl,gr,dmin,dmax,index_dmin_local,index_dmax_local,mynode,totalnodes);
    
    // we create a new vector that contains all the new fields. It is a support vector that will be swapped with the old one
    fields_vector new_fields_vect(N,vector<vector<vector<double>>>(dim_x,vector<vector<double>>(dim_y,vector<double>(dim_z))));

    for (int n=0; n <N; n++)
    {
        // Bulk processors evolving routine for k1
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)  // including the boundaries 
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    new_fields_vect[n][i][j][k] = fields_vect[n][i][j][k] + dt*(k1[n][i][j][k]+2*k2[n][i][j][k]+2*k3[n][i][j][k]+k4[n][i][j][k])/6.;
                //+artificial_diss(epsilon,ord,fields_vect,j,i,dx,dt));
                }
            }
        }
        //----communication(new_fields_vect,j,index_dmax_local,index_dmin_local,nitems,mynode,totalnodes,status);
        communication(new_fields_vect,Grid,n,index_dmax_local,index_dmin_local,mynode,totalnodes,status);
    }   
    

    // populate the ghost zone of the fields, we need them to calculate k1 in the next time step
    /*for (int j=0; j <N; j++)
    {
        ghost_point_extrapolation(new_fields_vect, t,dx,dt,j,gl,gr,dmin,dmax);
    }
    //cout<<"old "<<fields_vect.size()<<"new "<<new_fields_vect.size()<<endl;
    */
    ghost_point_extrapolation(new_fields_vect,Grid,S,t,dx,dt,gl,gr,dmin,dmax,index_dmin_local,index_dmax_local,mynode,totalnodes);
    return(new_fields_vect);
}



// initialize the grid function
grid initialization_grid(vector<double> dmin,vector<double> dmax,vector<double> dx,vector<int> gl, vector<int> gr,vector<int> ord,int dim)
{
    
    // setting the ghost point equal to the order of the artificial dissipator
    for(int d=0;d<dim;d++)
    {
        if (ord>gl)
        {
            gl[d] = ord[d];
        }
        if(ord>gr)
        {
            gr[d] = ord[d];
        }
    }
    // initializing the size of all the uni-dimensional grid
    
    vector<int> i_v(dim);
    for(int d=0;d<dim;d++)
    {
        if(dmax[d]!=dmin[d])
        {
            i_v[d] =  int( ((dmax[d]+dx[d]*double(gr[d]))-(dmin[d]-dx[d]*double(gl[d])) +dx[d]/2 ) / dx[d]) + 1;
        }
        else
        {
            i_v[d] = 1;
        }
            //cout<<i_v[d]<<endl;
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
                val1 = dmin[0]+dx[0]*(i-gl[0]);
                val2 = dmin[1]+dx[1]*(j-gl[1]);
                val3 = dmin[2]+dx[2]*(k-gl[2]);
                Grid[i][j].push_back({val1, val2, val3});
            }
        }
    }
    return(Grid);
} 

fields_vector initialize_fields(vector<double> dmin,vector<double> dmax,vector<double> dx,vector<int> gl, vector<int> gr,vector<int> ord,int dim,grid Grid,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> param_ic)
{
    int N = funcs.size();
    fields_vector new_fields;
    // setting the ghost point equal to the order of the artificial dissipator
    for(int d=0;d<dim;d++)
    {
        if (ord>gl)
        {
            gl = ord;
        }
        if(ord>gr)
        {
            gr = ord;
        }
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
                    new_fields[n][i][j].push_back(funcs[n](Grid[i][j][k],param_ic));
                }
            }
        }
    }
    return(new_fields);
} 

//fields_vect[n][i][j][k]
//,MPI_Status status, int totalnodes, int mynode,MPI_Request request
int print_f(fields_vector fields_vect, grid Grid, vector<double> dmin, vector<double> dmax, vector<double> dx, string name_file,string name_folder, vector<int> gl, vector<int> gr,int mynode, int totalnodes)
{
    int count_points = 0;
    
    int N = fields_vect.size();
    vector<int> S(dim);
    vector<double> index_dmin_local(dim),index_dmax_local(dim);
    for(int d=0;d<dim;d++)
    {
        //cout<<"max: "<<dmax[d]<<"min: "<<dmin[d]<<"dx :"<<dx[d]<<endl;
        if(dmax[d]!=dmin[d])
        {
            S[d] =  int( ((dmax[d]+dx[d]*double(gr[d]))-(dmin[d]-dx[d]*double(gl[d])) +dx[d]/2 ) / dx[d]) + 1;
        }
        else
        {
            S[d] = 1;
        }
        //cout<<"S"<<d<<": "<<S[d]<<endl;
        // defining the first and last spatial index of the domain, the actual process will work inside this range
        if (S[d]!=1)
        {
            index_dmin_local[d] = (gl[d]) + mynode * int((S[d]-gl[d]-gr[d])/totalnodes);
            //cout<<"mx from processor: "<<mynode<<" dmin local: dim "<<d<<" index "<<index_dmin_local[d]<<endl;
            if (mynode==totalnodes-1)
            {
                index_dmax_local[d] =  S[d]-gr[d];
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
            }
            else
            {
                index_dmax_local[d] = (gl[d]) + (mynode+1) * int((S[d]-gl[d]-gr[d])/totalnodes);
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
            }
        }
        else
        {
            index_dmin_local[d] = 0;
            index_dmax_local[d] = 1;
        }
    }
    
    
    // ! we parallelize on the first coordinate !
    for(int d=1;d<dim;d++)
    {
        if (S[d]!=1)
        {
            index_dmin_local[d] = (gl[d]) ;
            //cout<<"mx from processor: "<<mynode<<" dmin local: dim "<<d<<" index "<<index_dmin_local[d]<<endl;
           
            index_dmax_local[d] =  S[d]-gr[d];
                //cout<<"mx from processor: "<<mynode<<" dmax local: dim "<<d<<" index "<<index_dmax_local[d]<<endl;
        }
        else
        {
            index_dmin_local[d] = 0;
            index_dmax_local[d] = 1;
        }
    }
    
    
    // open the output file with name name_file without deleting its lines
    ofstream myfile_print_f;
    myfile_print_f.open (name_file,ios::app);
   
    
    // headers of the columns
    for(int d=0;d<dim;d++)
    {
        myfile_print_f<<"x_"<<to_string(d)<<",";
    }
    for (int i=0; i<fields_vect.size(); i++)
    {
        myfile_print_f << "field"<<to_string(i)<<",";
    }    
    myfile_print_f << "field"<<to_string(fields_vect.size()-1)<<"\n";
    
    
    // Each processor print only the spatial points that it is evolving
    
    // the processor zero print also the left ghost points
    if(mynode==0)
    {
        // we print the ghost point "before" the i=index_dmin_local plan[0]
        if(S[0]!=1)
        {
            //cout<<"we print the ghost point before the i=index_dmin_local plan[0]"<<endl;
            for(int i=0;i<gl[0];i++)
            {
                //cout<<"(gl)processor "<<mynode<<" index x_0: "<<i<<endl;
                for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
                {
                    //cout<<"(gl)processor "<<mynode<<" index x_1: "<<j<<endl;
                    for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                    {
                        count_points += 1;
                        //cout<<"(gl)processor "<<mynode<<" index x_2: "<<k<<endl;
                        for(int d=0;d<dim;d++)
                        {
                            myfile_print_f << Grid[i][j][k][d]<<",";
                        }
                        for(int n=0;n<N-1;n++)
                        {
                            myfile_print_f<< fields_vect[n][i][j][k]<<",";
                        }
                        
                        myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                        myfile_print_f<<"\n";
                    }
                }
            }
        }
    }
    // we print the ghost point "before" the j=index_dmin_local plan[1]
    if(S[1]!=1)
    {
        //cout<<"S[1] è diverso da 1"<<endl;
        //cout<<"we print the ghost point before the j=index_dmin_local plan[1]"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            //cout<<"(gl)processor "<<mynode<<" index x_0: "<<i<<endl;
            for(int j=0;j<gl[1];j++)
            {
                //cout<<"(gl)processor "<<mynode<<" index x_1: "<<j<<endl;
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    //cout<<"(gl)processor "<<mynode<<" index x_2: "<<k<<endl;
                    count_points += 1;
                    for(int d=0;d<dim;d++)
                    {
                        myfile_print_f << Grid[i][j][k][d]<<",";
                    }
                    for(int n=0;n<N-1;n++)
                    {
                        myfile_print_f << fields_vect[n][i][j][k]<<",";
                    }
                    myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                    myfile_print_f<<"\n";
                        
                }
            }
        }
    }
    // we print the ghost point "before" the k=index_dmin_local plan[2]
    if(S[2]!=1)
    {
        //cout<<"we print the ghost point before the k=index_dmin_local plan[2]"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=0;k<gl[2];k++)
                {
                    count_points += 1;
                    for(int d=0;d<dim;d++)
                    {
                        myfile_print_f << Grid[i][j][k][d]<<",";
                    }
                    for(int n=0;n<N-1;n++)
                    {
                        myfile_print_f << fields_vect[n][i][j][k]<<",";
                    }
                    myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                    myfile_print_f<<"\n";
                }
            }
        }
    }
    
    // Bulk processors printing routine
    for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
    {
        //cout<<"processor "<<mynode<<" index x_0: "<<i<<endl;
        for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
        {
            //cout<<"processor "<<mynode<<" index x_1: "<<j<<endl;
            for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
            {
                count_points += 1;
                //cout<<"processor "<<mynode<<" index x_2: "<<k<<endl;
                for(int d=0;d<dim;d++)
                {
                    myfile_print_f << Grid[i][j][k][d]<<",";
                }
                
                for(int n=0;n<N-1;n++)
                {
                    //cout<<"field index"<<n<<endl;
                    myfile_print_f << fields_vect[n][i][j][k]<<",";
                }
                myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                myfile_print_f<<"\n";
            }
        }
    }
                
    // RIGHT GHOST POINTS PRINTING
    // only the last processor print also the right ghost points
    if(mynode==totalnodes-1)
    {
        // we print the ghost point "after" the i=index_dmax_local[0] plan
        if(S[0]!=1)
        {
            //cout<<"we print the ghost point after the i=index_dmax_local[0] plan"<<endl;
            for(int i=index_dmax_local[0];i<S[0];i++)
            {
                //cout<<"last processor x_0 index: "<<i<<endl;
                for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
                {
                    //cout<<"last processor x_1 index: "<<j<<endl;
                    for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                    {
                        count_points += 1;
                        //cout<<"last processor x_2 index: "<<k<<endl;
                        for(int d=0;d<dim;d++)
                        {
                            myfile_print_f << Grid[i][j][k][d]<<",";
                        }
                        for(int n=0;n<N-1;n++)
                        {
                            myfile_print_f << fields_vect[n][i][j][k]<<",";
                        }
                        myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                        myfile_print_f<<"\n";
                    }
                }
            }
        }
    }
    // we print the ghost point "after" the j=index_dmax_local[1] plan
    if(S[1]!=1)
    {
        //cout<<"we print the ghost point after the j=index_dmax_local[1] plan"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmax_local[1];j<S[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    count_points += 1;
                    for(int d=0;d<dim;d++)
                    {
                        myfile_print_f << Grid[i][j][k][d]<<",";
                    }
                    for(int n=0;n<N-1;n++)
                    {
                        myfile_print_f << fields_vect[n][i][j][k]<<",";
                    }
                    myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                    myfile_print_f<<"\n";
                }
            }
        }
    }
    // we print the ghost point "after" the i=index_dmax_local plan[0]
    if(S[2]!=2)
    {
        //cout<<"we print the ghost point after the k=index_dmax_local[2] plan"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            //cout<<"last processor x_0 index: "<<i<<endl;
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                //cout<<"last processor x_1 index: "<<j<<endl;
                for(int k=index_dmax_local[2];k<S[2];k++)
                {
                    count_points += 1;
                    //cout<<"last processor x_2 index: "<<k<<endl;
                    for(int d=0;d<dim;d++)
                    {
                        myfile_print_f << Grid[i][j][k][d]<<",";
                    }
                    for(int n=0;n<N-1;n++)
                    {
                        myfile_print_f << fields_vect[n][i][j][k]<<",";
                    }
                    myfile_print_f << fields_vect[N-1][i][j][k]<<",";
                    myfile_print_f<<"\n";
                        
                }
            }
        }
    }
    
    
    myfile_print_f.close();
    return(count_points);
    
    

}

double we_PI_3D(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(Dx[0](fields_vec,0,i,j,k,dx[0])+Dx[1](fields_vec,1,i,j,k,dx[1])+Dx[2](fields_vec,2,i,j,k,dx[2]));
}

double we_PI_2D(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(Dx[0](fields_vec,0,i,j,k,dx[0])+Dx[1](fields_vec,1,i,j,k,dx[1]));
}

double we_PHIx(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(Dx[0](fields_vec,3,i,j,k,dx[0]));
}

double we_PHIy(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(Dx[1](fields_vec,3,i,j,k,dx[1]));
}

double we_PHIz(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(Dx[2](fields_vec,3,i,j,k,dx[2]));
}

double we_phi(int field_ind,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(fields_vec[3][i][j][k]);
}

double null_eq(int n,int i, int j, int k,fields_vector fields_vec, grid Grid,vector<double> dx,std::vector<double> param,double t,derivatives_vector Dx,double dt )
{
    return(0);
}

// derivative operator


double partial_dx(fields_vector fields_vec,int ind_field,int i, int j, int k,double dx)
{
    return((fields_vec[ind_field][i+1][j][k]
            -fields_vec[ind_field][i-1][j][k])  /2./dx);
}

double partial_dy(fields_vector fields_vec,int ind_field,int i, int j, int k,double dy)
{
    return((fields_vec[ind_field][i][j+1][k]
            -fields_vec[ind_field][i][j-1][k])  /2./dy);
}

double partial_dz(fields_vector fields_vec,int ind_field,int i, int j, int k,double dz)
{
    return((fields_vec[ind_field][i][j][k+1]
            -fields_vec[ind_field][i][j][k-1])  /2./dz);
}



// ---------------------- GHOST POINT EXTRAPOLATIONS ---------------------- //



//fields_vect[n][i][j][k]
void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid, vector<int> S, double t,vector<double> dx, double dt,vector<int> gl, vector<int> gr,vector<double> dmin,vector<double> dmax,vector<int> index_dmin_local, vector<int> index_dmax_local,int mynode,int totalnodes)
{   
    
    int N = field_vect.size();   
    
    // LEFT GHOST POINTS UPDATING //
    
    // only the processor zero update left ghost points of the axis 0
    if(mynode==0)
    {
        // we update the ghost point "before" the i=index_dmin_local plan[0]
        if(S[0]!=1)
        {
            for(int i=gl[0]-1;i>=0;i--)
            {
                for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
                {
                    for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                    {
                        for(int n=0;n<N;n++)
                        {
                            field_vect[n][i][j][k] = field_vect[n][i+1][j][k]-(-25./12.*field_vect[n][i+1][j][k]+4.*field_vect[n][i+2][j][k]-3.*field_vect[n][i+3][j][k]+4./3.*field_vect[n][i+4][j][k]-1./4.*field_vect[n][i+5][j][k])+1./2.*(2.*field_vect[n][i+1][j][k]-5.*field_vect[n][i+2][j][k]+4.*field_vect[n][i+3][j][k]-field_vect[n][i+4][j][k]);
                        }
                    }
                }
            }
        }
    }
    
    // we update the ghost point "before" the j=index_dmin_local plan[1]
    if(S[1]!=1)
    {
        //cout<<"we update the ghost point before the j=index_dmin_local plan[1]"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=gl[1]-1;j>=0;j--)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    for(int n=0;n<N;n++)
                    {
                        field_vect[n][i][j][k] = field_vect[n][i][j+1][k]-(-25./12.*field_vect[n][i][j+1][k]+4.*field_vect[n][i][j+2][k]-3.*field_vect[n][i][j+3][k]+4./3.*field_vect[n][i][j+4][k]-1./4.*field_vect[n][i][j+5][k])+1./2.*(2.*field_vect[n][i][j+1][k]-5.*field_vect[n][i][j+2][k]+4.*field_vect[n][i][j+3][k]-field_vect[n][i][j+4][k]);
                    }                        
                }
            }
        }
    }
    // we update the ghost point "before" the k=index_dmin_local plan[2]
    if(S[2]!=1)
    {
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=gl[2];k>=0;k--)
                {
                    for(int n=0;n<N;n++)
                    {
                        field_vect[n][i][j][k] = field_vect[n][i][j][k+1]-(-25./12.*field_vect[n][i][j][k+1]+4.*field_vect[n][i][j][k+2]-3.*field_vect[n][i][j][k+3]+4./3.*field_vect[n][i][j][k+4]-1./4.*field_vect[n][i][j][k+5])+1./2.*(2.*field_vect[n][i][j][k+1]-5.*field_vect[n][i][j][k+2]+4.*field_vect[n][i][j][k+3]-field_vect[n][i][j][k+4]);
                    }
                }
            }
        }
    }
    
    
    // RIGHT GHOST POINTS UPDATING //
    if(mynode==totalnodes-1)
    {
        // we print the ghost point "after" the i=index_dmax_local[0] plan
        if(S[0]!=1)
        {
            //cout<<"we print the ghost point after the i=index_dmax_local[0] plan"<<endl;
            for(int i=index_dmax_local[0];i<S[0];i++)
            {
                //cout<<"last processor x_0 index: "<<i<<endl;
                for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
                {
                    //cout<<"last processor x_1 index: "<<j<<endl;
                    for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                    {
                        //cout<<"last processor x_2 index: "<<k<<endl;
                        for(int n=0;n<N;n++)
                        {
                            field_vect[n][i][j][k] = field_vect[n][i-1][j][k] +(+25./12.*field_vect[n][i-1][j][k] -4.*field_vect[n][i-2][j][k] +3.*field_vect[n][i-3][j][k] -4./3.*field_vect[n][i-4][j][k] +1./4.*field_vect[n][i-5][j][k] +1./2.*(2.*field_vect[n][i-1][j][k] -5.*field_vect[n][i-2][j][k] +4.*field_vect[n][i-3][j][k] -field_vect[n][i-4][j][k] ));
                        }
                    }
                }
            }
        }
    }
    // we print the ghost point "after" the j=index_dmax_local[1] plan
    if(S[1]!=1)
    {
        //cout<<"we print the ghost point after the j=index_dmax_local[1] plan"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmax_local[1];j<S[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    for(int n=0;n<N;n++)
                    {
                        field_vect[n][i][j][k] = field_vect[n][i][j-1][k] +(+25./12.*field_vect[n][i][j-1][k] -4.*field_vect[n][i][j-2][k] +3.*field_vect[n][i][j-3][k] -4./3.*field_vect[n][i][j-4][k] +1./4.*field_vect[n][i][j-5][k] +1./2.*(2.*field_vect[n][i][j-1][k] -5.*field_vect[n][i][j-2][k] +4.*field_vect[n][i][j-3][k] -field_vect[n][i][j-4][k] ));
                    }
                }
            }
        }
    }
    // we print the ghost point "after" the i=index_dmax_local plan[0]
    if(S[2]!=2)
    {
        //cout<<"we print the ghost point after the k=index_dmax_local[2] plan"<<endl;
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            //cout<<"last processor x_0 index: "<<i<<endl;
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                //cout<<"last processor x_1 index: "<<j<<endl;
                for(int k=index_dmax_local[2];k<S[2];k++)
                {
                    for(int n=0;n<N;n++)
                    {
                        field_vect[n][i][j][k] = field_vect[n][i][j][k-1] +(+25./12.*field_vect[n][i][j][k-1] -4.*field_vect[n][i][j][k-2] +3.*field_vect[n][i][j][k-3] -4./3.*field_vect[n][i][j][k-4] +1./4.*field_vect[n][i][j][k-5] +1./2.*(2.*field_vect[n][i][j][k-1] -5.*field_vect[n][i][j][k-2] +4.*field_vect[n][i][j][k-3] -field_vect[n][i][j][k-4] ));
                    }
                }
            }
        }
    }
   
}

void communication(fields_vector &fields_vect,grid Grid,int n,vector<int> index_dmax_local,vector<int> index_dmin_local,int mynode,int totalnodes,MPI_Status status)
{
    int dim_x = Grid.size(),dim_y=Grid[0].size(),dim_z=Grid[0][0].size();
   
    int nitems = 2*dim_y*dim_z;
    double message_out[2][dim_y][dim_z];
    double message_in[2][dim_y][dim_z];
    double message_out_left[2][dim_y][dim_z];
    double message_out_right[2][dim_y][dim_z];
    double message_in_left[2][dim_y][dim_z];
    double message_in_right[2][dim_y][dim_z];
    
    // note: we are parallelizing in the first dimension
    // note: the node zero will communicate only it's right surface
    // NODE ZERO COMMUNICATION //
    if(mynode==0)
    {
        // filling the out put message
        for(int i=index_dmax_local[0]-2;i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    message_out[i-(index_dmax_local[0]-2)][j][k] = fields_vect[n][i][j][k] ;
                }
            }
        }
        // communication between node zero and node 1
        MPI_Sendrecv(&message_out,nitems, MPI_DOUBLE, mynode+1,
        0, &message_in, nitems, MPI_DOUBLE, mynode+1,  0, MPI_COMM_WORLD, &status);
        // extraction of the input message
        for(int i=index_dmax_local[0];i<index_dmax_local[0]+2;i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    fields_vect[n][i][j][k] = message_in[i-(index_dmax_local[0])][j][k];
                }
            }
            
        }
    }
    // LAST NODE COMMUNICATION //
    else if(mynode==totalnodes-1)
    {
        // filling the out put message
        for(int i=index_dmin_local[0];i<index_dmin_local[0]+2;i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    message_out[i-(index_dmin_local[0])][j][k] = fields_vect[n][i][j][k] ;
                }
            }
        }
        MPI_Sendrecv(&message_out,nitems, MPI_DOUBLE, mynode-1,
        0, &message_in, nitems, MPI_DOUBLE, mynode-1,  0, MPI_COMM_WORLD, &status);
        // extraction of the input message
        for(int i=index_dmin_local[0]-2;i<index_dmin_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    fields_vect[n][i][j][k] = message_in[i-(index_dmin_local[0]-2)][j][k];
                }
            }
            
        }
    }
    // BULK NODE COMMUNICATION //
    else
    {
        // BUILDING OUTPUT MESSAGES 
        // FILLING the LEFT out put message
        for(int i=index_dmin_local[0];i<index_dmin_local[0]+2;i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    message_out_left[i-(index_dmin_local[0])][j][k] = fields_vect[n][i][j][k] ;
                }
            }
        }
        // FILLING the RIGHT out put message
        for(int i=index_dmax_local[0]-2;i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    message_out_right[i-(index_dmax_local[0]-2)][j][k] = fields_vect[n][i][j][k] ;
                }
            }
        }
        // SENDING THE MESSAGES
        MPI_Send(&message_out_left,nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD);
        MPI_Send(&message_out_right,nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD);
        
        //RECEIVING THE MESSAGES
        MPI_Recv(&message_in_left,nitems,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD, &status);
        MPI_Recv(&message_in_right,nitems,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD, &status);
        // EXTRACTION of the LEFT input message
        for(int i=index_dmin_local[0]-2;i<index_dmin_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    fields_vect[n][i][j][k] = message_in_left[i-(index_dmin_local[0]-2)][j][k];
                }
            }
            
        }
        // EXTRACTION of the RIGHT input message
        for(int i=index_dmax_local[0];i<index_dmax_local[0]+2;i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for( int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
                    fields_vect[n][i][j][k] = message_in_right[i-(index_dmax_local[0])][j][k];
                }
            }
            
        }
    }
}


double initial_gauss(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( init_param[0]*exp(-pow((x[0]*dev_std),2)-pow((x[1]*dev_std),2)-pow((x[2]*dev_std),2)) );
}

double initial_gauss_PHIx(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( -2*x[0]*init_param[0]*exp(-pow((x[0]*dev_std),2)-pow((x[1]*dev_std),2)-pow((x[2]*dev_std),2)) );
}

double initial_gauss_PHIy(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( -2*x[1]*init_param[0]*exp(-pow((x[0]*dev_std),2)-pow((x[1]*dev_std),2)-pow((x[2]*dev_std),2)) );
}

double initial_null(vector<double> x,vector<double> init_param)
{
    return(0);
}
