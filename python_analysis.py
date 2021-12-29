# support function to perform data analysis #

import numpy as np
import matplotlib.pyplot as plt

    


def gauss_PI(x,t):
    return(-2*(x+t)*np.exp(-pow(x+t,2))+2*(x-t)*np.exp(-pow(x-t,2)) )
    
def gauss_phi(x,t):
    return(np.exp(-pow(x+t,2))+np.exp(-pow(x-t,2)))
    
def sin_func(x,t):
    return(np.sin(x+t))

def line(x,t):
    return(4*t)

def conv_test(vect1,vect2,theo,t):
    field1 = tuple(vect1['field0'])
    field2 = tuple(vect2['field0'])
    
    dx1 = vect1['x'][1]-vect1['x'][0]
    dx2 = vect2['x'][1]-vect2['x'][0]
    dmin = vect2['x'][0]
    dmax = vect1['x'][vect1['x'].size-1]
    
    theo_vect1 = initialize_func_vect(theo,vect1['x'] ,dx1,t)
    #plt.plot(field1)
    #plt.plot(theo_vect1)
    norm_diff_1 = norm(np.subtract(field1,theo_vect1),dx1)
    #plt.plot(vect1['x'],np.subtract(field1,theo_vect1))
    
    theo_vect2 = initialize_func_vect(theo,vect2['x'] ,dx2,t)
    norm_diff_2 = norm(np.subtract(field2,theo_vect2),dx2)
    #plt.plot(vect2['x'],np.subtract(field2,theo_vect2))
    return(np.log(norm_diff_1/norm_diff_2)/np.log(2))
    
def self_conv_test(vect1, vect2, vect3,gl,gr):
    field1 = tuple(vect1['field0'][gl:-gr])
    field2 = tuple(vect2['field0'][gl:-gr])
    field3 = tuple(vect3['field0'][gl:-gr])
    x1 = tuple(vect1['x'][gl:-gr])
    x2 = tuple(vect2['x'][gl:-gr])
    dx1 = vect1['x'][1]-vect1['x'][0]
    dx2 = vect2['x'][1]-vect2['x'][0]
    #dx3 = vect3['x'][1]-vect3['x'][0]
    dmin = vect1['x'][0]
    dmax = vect1['x'][vect1['x'].size-1]
    norm_diff_1 = norm(np.subtract(field1,field2[::2]),x1,dx1)
    norm_diff_2 = norm(np.subtract(field2,field3[::2]),x2,dx2)
    return(np.log(norm_diff_1/norm_diff_2)/np.log(2.0))

def self_conv_test_spherical(vect1, vect2, vect3,gl,gr,field):
    field1 = tuple(vect1[field][gl:-gr])
    field2 = tuple(vect2[field][gl:-gr])
    field3 = tuple(vect3[field][gl:-gr])
   
    
    dx1 = vect1['x'][1]-vect1['x'][0]
    dx2 = vect2['x'][1]-vect2['x'][0]
    x1 = tuple(vect1['x'][gl:-gr])
    x2 = tuple(vect2['x'][gl:-gr])
    #dx3 = vect3['x'][1]-vect3['x'][0]
    norm_diff_1 = spherical_norm(np.subtract(field1,field2[::2]),x1,dx1)
    norm_diff_2 = spherical_norm(np.subtract(field2,field3[::2]),x2,dx2)
    return(np.log(norm_diff_1/norm_diff_2)/np.log(2.0))


def self_conv_test_pw(vect1, vect2, vect3,gl,gr,field):
    field1 = tuple(vect1[field][gl:-gr])
    field2 = tuple(vect2[field][gl:-gr])
    field3 = tuple(vect3[field][gl:-gr])
    diff_1 = np.subtract(field1,field2[::2])
    diff_2 = np.subtract(field2[::2],field3[::4])
    return(diff_1,4*diff_2)

def norm(vect,x,dx):
    norm = 0
    for i in range(len(vect)):
        if (vect[i] == 'Nan'):
            return('Nan')
        norm = norm + vect[i]**2 *dx
    norm = norm**(0.5)
    return (norm)

def spherical_norm(vect,x,dx):
    norm = 0
    for i in range(len(vect)):
        if (vect[i] == 'Nan'):
            return('Nan')
        norm = norm + vect[i]**2 *x[i]**2 *dx
    norm = norm**(0.5)
    return (norm)
    
    
def sin_func(x,t):
    return(np.sin(x+t))

def parab(x,t):
    return((x+t)**2)

def cube(x,t):
    return((x+t)**3)

def iperb(x,t):
    return(1/(x+t+1))

def spherical_we_gaussian_solution(x,t,a):
    return(-a/4*(np.exp(-1*(t+x)**2)-np.exp(-1*(t-x)**2) )  /x)

def model1_gaussian_solution(x,t,a):
    return(np.log(1+(-a/4*(np.exp(-1*(t+x)**2)-np.exp(-1*(t-x)**2) )  /x) ))


def model3_gaussian_solution(x,t,a):
    return(np.sin(np.log(1+(-a/4*(np.exp(-1*(t+x)**2)-np.exp(-1*(t-x)**2) )  /x) )))
    
def initialize_func_vect(func,domain,dx,t):
    vect = []
    for x in domain:
        vect.append(func(x,t))
    return(vect)

# big_DF[ind_run][ind_dx][ind_time][ind_field]
def runs_maximums_vector(vectors,ind_dx,field):
    maximums_vector_run = []
    for i in range(len(vectors)): # cycle on the runs
        maximums_vector = []
        for j in range(len(vectors[i][ind_dx])): # cycle on the times
            value = max(np.absolute(vectors[i][ind_dx][j][field]))
            maximums_vector.append(value)
        maximums_vector_run.append(max(maximums_vector))
    return(maximums_vector_run)

def amplitudes_vector_creator(minimum,maximum,critical,dx):
    xrange = np.arange(minimum,maximum,dx)
    yrange = critical-np.exp(-xrange)
    return(yrange)

def self_similar_cordinates(space_vector, time_vector):
    for t in range(len(time_vector)):
        T.append(-np.log(time_vector[t]))
        x = []
        for i in range(len(big_DF[0][2][0]['x'])):
            x.append(space_vector[i]/time_vector[t])
        X.append(x)
    return(X,T)
### script to call in the main ###
"""

# create a vector of the convergence test. 
# conv_test_vect = []
for i,t in enumerate(times):
    conv_test_vect.append(conv_test(DF2[i],DF3[i],gauss_PI,t))
#print(conv_test_vect)
plt.plot(times,conv_test_vect)
#plt.xlim(40,80)
#plt.ylim(1.8,2.2)
#plt.vlines(47,2,0,alpha=0.2)
plt.grid()
"""