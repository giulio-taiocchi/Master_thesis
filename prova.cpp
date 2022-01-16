// Epa! //
#include <iostream>
#include <fstream>
#include<vector>
#include <math.h>
#include <chrono>
#include "boost/multi_array.hpp"
using namespace std;

typedef vector<vector<vector< vector<double> >>> grid;
typedef vector<vector<vector< vector<double> >>> fields_vector;

grid initialization_grid(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim);

grid initialize_fields(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim,grid Grid,std::vector<double(*)(vector<double>,vector<double>)> funcs,vector<double> &param_ic);

double initial_gauss(vector<double> x,vector<double> init_param);

print_f(fields_vector fields_vect, grid Grid vector<double> &dmin, vector<double> &dmax, vector<double> dx, string name_file,string name_folder, int gl, int gr);

void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax);

int main() 
{
    int mynode, totalnodes;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
    
    int dim = 3;
    vector<double> dmin {-1.,-1.,0};
    vector<double> dmax {1.,1.,0};
    vector<double> dx {0.5,0.5,0.};
    int gl = 2, gr = 2,ord = 2;
    int step_to_save = 1;
    double integration_interval = 1;
    
    vector<double> ic_par = {1};
    
    
    string file_path = "./data/prova/";
    string name_file = file_path + "processor_"+to_string(mynode)+"_ampl_"+to_string(ic_par[0])+"_eps"+to_string(0)+"_dx_"+to_string(dx[0])+"steps"+to_string(step_to_save)+"last_time"+to_string(integration_interval)+".csv";
    
    vector<double(*)(vector<double>, vector<double>)> initial_conditions;
    initial_conditions.push_back(&initial_gauss);    //PI1
    initial_conditions.push_back(&initial_gauss);        //PHI1       
    initial_conditions.push_back(&initial_gauss);        //phi1
    
    grid Grid= initialization_grid(dmin,dmax,dx, gl,gr,ord,dim);
    
    fields_vector new_fields_vector = initialize_fields(dmin,dmax,dx, gl,  gr, ord, dim, Grid,initial_conditions,ic_par);
    
    for(int n=0;n<initial_conditions.size();n++)
    {
        for(int i=0;i<Grid.size();i++)
        {
            for(int j=0;j<Grid[i].size();j++)
            {
                for(int k=0;k<Grid[i][j].size();k++)
                {
                    cout<<"point: ";
                    for(int d=0;d<dim;d++)
                    {
                        cout<<Grid[i][j][k][d]<<"  ";
                    }
                    cout<<"\nfields: ";
                    cout<<new_fields_vector[n][i][j][k];
                    cout<<endl;
                }
                cout<<endl;
            }
            cout<<endl;
        }
    }
    print_f(fields_vect, Grid dmin, dmax, dx, name_file, name_file,  gl,  gr)
    return 0;
} 



// initialize the grid function
grid initialization_grid(vector<double> &dmin,vector<double> &dmax,vector<double> &dx,int gl, int gr,int ord,int dim)
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

//fields_vect[n][i][j][k]
//,MPI_Status status, int totalnodes, int mynode,MPI_Request request
void print_f(fields_vector fields_vect, grid Grid vector<double> &dmin, vector<double> &dmax, vector<double> dx, string name_file,string name_folder, int gl, int gr)
{
    int dim = 3;
    int N = fields_vect.size()
    vector<int> S(dim);
    vector<double> index_dmin_local(dim),index_dmax_local(dim);
    for(int d=0;d<dim;d++)
    {
        S[d] =  int( ((dmax[d]+dx[d]*double(gr))-(dmin[d]-dx[d]*double(gl)) +dx[d]/2 ) / dx[d]) + 1;
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
        for(int i=0;i<gl+1;i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
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
        // we print the ghost point "before" the j=index_dmin_local plan[1]
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=0;j<gl+1;j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
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
        // we print the ghost point "before" the k=index_dmin_local plan[2]
        for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=0;k<gl+1;k++)
                {
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
    
    for(int i=index_dmin_local[0];i<index_dmax_local[0];i++)
    {
        for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
        {
            for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
            {
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
                
                
    // the last processor print also the right ghost points
    if(mynode==totalnodes-1)
    {
        // we print the ghost point "after" the i=index_dmax_local[0] plan
        for(int i=index_dmax_local[0];j<S[0].size();j++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
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
        // we print the ghost point "after" the j=index_dmax_local[1] plan
        for(int i=index_dmin_local[0];j<index_dmax_local[0].size();j++)
        {
            for(int j=index_dmax_local[1];j<S[1];j++)
            {
                for(int k=index_dmin_local[2];k<index_dmax_local[2];k++)
                {
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
        // we print the ghost point "after" the i=index_dmax_local plan[0]
        for(int i=index_dmin_local[0];j<index_dmax_local[0].size();j++)
        {
            for(int j=index_dmin_local[1];j<index_dmax_local[1];j++)
            {
                for(int k=index_dmax_local[2];k<S[2];k++)
                {
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
    myfile_print_f.close();
}

void ghost_point_extrapolation_4_ord(fields_vector &field_vect,grid Grid,double t,vector<double> dx, double dt, int j,int gl, int gr,vector<double> &dmin,vector<double> &dmax);

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

double initial_gauss(vector<double> x,vector<double> init_param)
{
    double dev_std = 1.;
    return( init_param[0]*exp(-pow((x[0]*dev_std),2)-pow((x[1]*dev_std),2)-pow((x[2]*dev_std),2)) );
}
