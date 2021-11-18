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
   
    
    dx1 = vect1['x'][1]-vect1['x'][0]
    dx2 = vect2['x'][1]-vect2['x'][0]
    #dx3 = vect3['x'][1]-vect3['x'][0]
    dmin = vect1['x'][0]
    dmax = vect1['x'][vect1['x'].size-1]
    norm_diff_1 = norm(np.subtract(field1,field2[::2]),dx1)
    norm_diff_2 = norm(np.subtract(field2,field3[::2]),dx2)
    return(np.log(norm_diff_1/norm_diff_2)/np.log(2.0))

def self_conv_test_pw(vect1, vect2, vect3,gl,gr):
    field1 = tuple(vect1['field0'][gl:-gr])
    field2 = tuple(vect2['field0'][gl:-gr])
    field3 = tuple(vect3['field0'][gl:-gr])
    dx1 = vect1['x'][1]-vect1['x'][0]
    dx2 = vect2['x'][1]-vect2['x'][0]
    #dx3 = vect3['x'][1]-vect3['x'][0]
    dmin = vect1['x'][0]
    dmax = vect1['x'][vect1['x'].size-1]
    diff_1 = np.subtract(field1,field2[::2])
    diff_2 = np.subtract(field2[::2],field3[::4])
    return(diff_1,4*diff_2)

def norm(vect,dx):
    norm = 0
    for i in range(len(vect)):
        norm = norm + vect[i]**2 *dx
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

def initialize_func_vect(func,domain,dx,t):
    vect = []
    for x in domain:
        vect.append(func(x,t))
    return(vect)