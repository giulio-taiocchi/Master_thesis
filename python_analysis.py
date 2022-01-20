# support function to perform data analysis #

import numpy as np
import pandas as pd
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

#def model3_PI_gaussian_solution(x,t,a):
 #   return(np.sin(np.log(1+(-a/4*(np.exp(-1*(t+x)**2)-np.exp(-1*(t-x)**2) )  /x) )))
    
    
def hyperbolic_chi_we_solution(r,t,a,ds,s):
    return(-a*np.exp(-ds**2*(-r+t)**2)+a*np.exp(-ds**2*(r**3+r*s**2-r**2*t+s**2*t)**2/(r**2-s**2)**2)*(1+r**2/(1-r**2/s**2)**2)**0.5*(1-r**2/s**2)/r)
    
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

def read_parallel_data(d_max,d_min,gl,gr,domain_lenght,h1,h2,h3,number_of_proc,number_steps,names):
    big_DF_sup=[]
    S = []
    for dx in [h1,h2,h3]:
        S.append(int( ((d_max+dx*gr)-(d_min-dx*gl) +dx/2 ) / dx) + 1)
    print(S)
    N_point1 = int((S[0]-2-gl-gr)/number_of_proc)
    N_point2 = int((S[1]-2-gl-gr)/number_of_proc)
    N_point3 = int((S[2]-2-gl-gr)/number_of_proc)
    # the processor zero has gl+1 ponts more: the left ghost points and the left boundary
    N_point1_zero = N_point1+gl+1
    N_point2_zero = N_point2+gl+1
    N_point3_zero = N_point3+gl+1
    
    N_point1_last = S[0]-1-gr-gl-1-(number_of_proc-1)*int((S[0]-2-gl-gr)/number_of_proc)+gr+1
    N_point2_last = S[1]-1-gr-gl-1-(number_of_proc-1)*int((S[1]-2-gl-gr)/number_of_proc)+gr+1
    N_point3_last = S[2]-1-gr-gl-1-(number_of_proc-1)*int((S[2]-2-gl-gr)/number_of_proc)+gr+1
    print(N_point1,N_point1_zero,N_point1_last)
    '''
    N_point1 = int((int(domain_lenght/h1-2))/number_of_proc)
    N_point2 = int((int(domain_lenght/h2-2))/number_of_proc)
    N_point3 = int((int(domain_lenght/h3-2))/number_of_proc)
    
    N_point1_last = N_point1+6
    N_point2_last = N_point2+6
    N_point3_last = N_point3+6
    '''
    total_point_zero = (N_point1_zero+1)*(number_steps+1)+ (N_point2_zero+1)*(number_steps+1)+(N_point3_zero+1)*number_steps+N_point3_zero
    total_point_others = (N_point1+1)*(number_steps+1)+ (N_point2+1)*(number_steps+1)+(N_point3+1)*number_steps+N_point3
    total_point_last = (N_point1_last+1)*(number_steps+1)+ (N_point2_last+1)*(number_steps+1)+(N_point3_last+1)*number_steps+N_point3_last
    total_points_vec = [total_point_zero]
    for n_proc in range (1,number_of_proc-1):
        total_points_vec.append(total_point_last)
    total_points_vec.append(total_point_last)
    DF_concatenated = []
    for n in range (0, len(names[:]),number_of_proc):
        DF_support = pd.read_csv(names[n] )
        #print(DF_support['x'][0]*2)
        for n_proc in range (1,number_of_proc):
            DF_support = pd.concat([DF_support,pd.read_csv(names[n+n_proc])])
        DF_concatenated.append(DF_support)

    for n in range (0, len(DF_concatenated[:])):
        for n_proc in range (0,number_of_proc):
            DF1, DF2, DF3 = [], [], []
            for i in range (0,number_steps):
                if(n_proc==0):
                    DF1.append( DF_concatenated[n][(N_point1_zero+1)*i:(N_point1_zero+1)*i+N_point1_zero].astype(float))
                    DF2.append( DF_concatenated[n][(N_point1_zero+1)*(number_steps+1)+ (N_point2_zero+1)*i : (N_point1_zero+1)*(number_steps+1)+ (N_point2_zero+1)*i+N_point2_zero].astype(float) )
                    DF3.append( DF_concatenated[n][(N_point1_zero+1)*(number_steps+1)+ (N_point2_zero+1)*(number_steps+1)+(N_point3_zero+1)*i  :  (N_point1_zero+1)*(number_steps+1)+ (N_point2_zero+1)*(number_steps+1)+(N_point3_zero+1)*i+N_point3_zero].astype(float))
                if(n_proc==number_of_proc-1):
                    DF1.append(DF_concatenated[n][(N_point1_last+1)*i+total_point_others*(n_proc-1)+total_point_zero:(N_point1_last+1)*i+N_point1_last+total_point_others*(n_proc-1)+total_point_zero].astype(float))
                    DF2.append( DF_concatenated[n][(N_point1_last+1)*(number_steps+1)+ (N_point2_last+1)*i +total_point_others*(n_proc-1)+total_point_zero: (N_point1_last+1)*(number_steps+1)+ (N_point2_last+1)*i+N_point2_last+total_point_others*(n_proc-1)+total_point_zero].astype(float) )
                    DF3.append( DF_concatenated[n][(N_point1_last+1)*(number_steps+1)+ (N_point2_last+1)*(number_steps+1)+(N_point3_last+1)*i+total_point_others*(n_proc-1)+total_point_zero:(N_point1_last+1)*(number_steps+1)+ (N_point2_last+1)*(number_steps+1)+(N_point3_last+1)*i+N_point3_last+total_point_others*(n_proc-1)+total_point_zero].astype(float))
                if(n_proc==1 and n_proc!=number_of_proc-1):
                    DF1.append(DF_concatenated[n][(N_point1+1)*i+total_point_zero:(N_point1+1)*i+N_point1+total_point_zero].astype(float))
                    DF2.append(DF_concatenated[n][(N_point1+1)*(number_steps+1)+ (N_point2+1)*i +total_point_zero: (N_point1+1)*(number_steps+1)+ (N_point2+1)*i+N_point2+total_point_zero].astype(float) )
                    DF3.append(DF_concatenated[n][(N_point1+1)*(number_steps+1)+ (N_point2+1)*(number_steps+1)+(N_point3+1)*i  +total_point_zero:(N_point1+1)*(number_steps+1)+ (N_point2+1)*(number_steps+1)+(N_point3+1)*i+N_point3+total_point_zero].astype(float))
                if(n_proc!=0 and n_proc!=number_of_proc-1 and n_proc !=1):
                    DF1.append(DF_concatenated[n][(N_point1+1)*i+total_point_others*(n_proc-1)+total_point_zero:(N_point1+1)*i+N_point1+total_point_others*(n_proc-1)+total_point_zero].astype(float))
                    DF2.append(DF_concatenated[n][(N_point1+1)*(number_steps+1)+ (N_point2+1)*i+total_point_others*(n_proc-1)+total_point_zero : (N_point1+1)*(number_steps+1)+ (N_point2+1)*i+N_point2+total_point_others*(n_proc-1)+total_point_zero].astype(float).astype(float) )
                    DF3.append(DF_concatenated[n][(N_point1+1)*(number_steps+1)+ (N_point2+1)*(number_steps+1)+(N_point3+1)*i +total_point_others*(n_proc-1)+total_point_zero :  (N_point1+1)*(number_steps+1)+ (N_point2+1)*(number_steps+1)+(N_point3+1)*i+N_point3+total_point_others*(n_proc-1)+total_point_zero].astype(float))
            big_DF_sup.append([DF1,DF2,DF3])
        print("run:"+str(n)+" ->"+names[n*number_of_proc]+ " added")
    print("number of runs:",len(big_DF_sup))
    new_big_DF = []
    for n in range (0, len(names[:]),number_of_proc):
        resolution_vector = []
        for i in range(0,3):
            time_steps_vector=[]
            for t in range(0,number_steps):
                time_steps_vector.append(pd.concat([big_DF_sup[n][i][t],big_DF_sup[n+1][i][t],big_DF_sup[n+2][i][t],big_DF_sup[n+3][i][t]]).reset_index())
            resolution_vector.append(time_steps_vector)
        new_big_DF.append(resolution_vector)
    return(new_big_DF)

def read_3D_parallel_data(names,number_of_proc):
    Big_Fields = []
    Big_Grid = []
    for i in range(0,len(names[:]),number_of_proc):
        Grid = []
        Fields = []
        fields = ["field0","field1","field2"]
        N_points = []
        DF_support = []
        for n in range (0, number_of_proc):
            DF_support.append( pd.read_csv(names[i+n]))
            N_points.append(int(np.array(DF_support[n]["x_0"])[-1]))
        Grid = np.array(DF_support[0][["x_0","x_1","x_2"]])[0:N_points[0]].astype(float)
        for n in range (1, number_of_proc):
            Grid_supp = np.array(DF_support[n][["x_0","x_1","x_2"]])[0:N_points[n]].astype(float)
            Grid = np.append(Grid,Grid_supp,axis=0)
        print("Grid shape ",Grid.shape)
        Fields = np.array(DF_support[0][fields])[0:N_points[0]].astype(float)
        for n in range (1, number_of_proc):
            Fields_supp = np.array(DF_support[n][fields])[0:N_points[n]].astype(float)
            Fields = np.append(Fields,Fields_supp,axis=0)
        Fields.shape
        Big_Grid.append(Grid)
        Big_Fields.append(Fields)
    return(Big_Grid)

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