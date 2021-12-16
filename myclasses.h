#include <math.h>
#include <fstream>
#include <sstream>
#include "spline.h"
using namespace std;


class layer
{
    // Access specifier
    public:
 
    // Data Members
    double dmin, dmax,dx;
    vector<double> grid;
    double time,dt;
    int n_field
    vector<vector<double>> grid_functions;
    layer(double dmin_in, double dmax_in, double dx_in, double dt_in, int n_field_in,vector<vector<double>> grid_functions_in)
    {
        dmin = dmin_in;
        dmax = dmax_in;
        dx = dx_in;
        dt = dt_in;
        n_field = n_field_in;
        vector<vector<double>> grid_functions = grid_functions_in;
    }
    
 
    // Member Functions()
    void printname()
    {
       cout << "Geekname is: " << geekname;
    }
}; 
