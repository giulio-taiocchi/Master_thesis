# support function to perform data analysis #
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from python_analysis import *
import matplotlib.animation as animation
from numpy import random
from sklearn.linear_model import LinearRegression
from scipy.signal import argrelextrema


    


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
    
    
#def hyperbolic_chi_we_solution(r,t,a,ds,s):
    #return(-a*np.exp(-ds**2*(-r+t)**2)+a*np.exp(-ds**2*(r**3+r*s**2-r**2*t+s**2*t)**2/(r**2-s**2)**2)*(1+r**2/(1-r**2/s**2)**2)**0.5*(1-r**2/s**2)/r) 
 

def hyperbolic_chi_we_solution(r,t,a,ds,s):
    return((a*(r-t)*np.exp(-1*ds**2*(r-t)**2)+a*(-r+t+2*r/(1-r**2/s**2))*np.exp(-1*ds**2*(-r+t+2*r/(1-r**2/s**2))**2 ))*(1+r**2/(1-r**2/s**2)**2)**0.5*(1-r**2/s**2)/2/r)

def spherical_we_solution(r,t,a,ds):
    return( (a*(r-t)*np.exp(-(ds*(r-t))**2)+a*(r+t)*np.exp(-(ds*(r+t))**2))/2/r )

def initial_test(r,t,a,ds,s):
    return(a*(1+r**2/(1-r**2/s**2)**2)**0.5*np.exp(-ds**2*(r/(1-r**2/s**2))**2))

def initialize_func_vect(func,domain,dx,t):
    vect = []
    for x in domain:
        vect.append(func(x,t))
    return(vect)

# MODEL 1 WITH CHARACTERISTIC VARIABLES, hyperboloidal compactification, Chi rescaling, THEORETICAL SOLUTIONS 

def hyperbolic_chi_m1_charvar_solution_Psi(r,t,a,B,ds,s):
    A = a
    R = r/(1-r**2/s**2)
    T = t+R-r
    Chi = (1+R**2)**0.5
    return( np.log(1+B/(2*R)*(A*np.exp(-ds**2*(R-T)**2)*(R-T)+A*np.exp(-ds**2*(R+T)**2)*(R+T)))*Chi )


def hyperbolic_chi_m1_charvar_solution_PsiPlus(r,t,a,B,ds,s):
    A = a
    R = r/(1-r**2/s**2)
    T = t+R-r
    Chi = (1+R**2)**0.5
    J = R-T
    K = np.exp(4*ds**2*R*T)
    L = (R+T)**2
    return(  (A*B*(J-J*K-4*ds**2*R*L)*Chi**2)/(R*(2*np.exp(ds**2*L)*R+A*B*(J*K+R+T)))  )

def hyperbolic_chi_m1_charvar_solution_PsiMinus(r,t,a,B,ds,s):
    A = a
    R = r/(1-r**2/s**2)
    T = t+R-r
    Chi = (1+R**2)**0.5
    J = R-T
    K = np.exp(4*ds**2*R*T)
    L = (R+T)
    return(  A*B*(L+K*((-1+2*ds*J)*(1+2*ds*J)*R-T))*Chi )/(R*(A*B*(J*K+L)+2*np.exp(ds**2*L**2)*R))


def hyperbolic_chi_m3_charvar_solution_Psi2Plus(r,t,a,B,C,ds,s):
    A = a
    R = r/(1-r**2/s**2)
    T = t+R-r
    Chi = (1+R**2)**0.5
    J = R-T
    K = np.exp(4*ds**2*R*T)
    L = (R+T)
    M = R**2+T**2
    N = B*A*(np.exp(-ds**2*J**2)*J+np.exp(-ds**2*L**2)*L)/2/R
    P = np.log(1+N)/C
    Q = np.log(1+(A*np.exp(-ds**2*J**2)*J+A*np.exp(-ds**2*L**2)*L)/2/R)/C
    return(A*np.exp(ds**2*J**2)*np.sin(Q)*(J*K+(-1+4*ds**2*L**2)*R+T)*Chi**2 )/R/(A*np.exp(ds**2*J**2)*(J*K+L)+2*np.exp(2*ds**2*M)*R)

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

def names_generator(data,epsilon,amplitude_vector,number_of_proc,model,dx,number_steps,range_time):
    names = []
    for i in data:
        for e in range(len(epsilon)):
            for a in range(len(amplitude_vector)):
                for n in range(number_of_proc):
                    names.append("./data"+str(i)+"/processor_"+                                                                                                str(n)+"_ampl_"+str(format(amplitude_vector[a],'.6f'))+"_eps"+str(format(epsilon[e],'.6f'))+"_dx_"+str(format(dx,'.6f'))+"steps"+str(number_steps)+"last_time"+str(format(range_time, '.6f'))+".csv") 
    return(names)
    
    
### PLOT FUNCTIONS ###
def plot_pw_convergence(big_DF,gl,gr,self_conv_test_vect_pw,model,data,field,epsilon,n_ind,dx,number_steps,ylim_inf,ylim_sup):
    fig, ax = plt.subplots()

    line1, = ax.plot(big_DF[n_ind][0][1]['x'][gl:-gr], self_conv_test_vect_pw[1][0],'-.')
    line2, = ax.plot(big_DF[n_ind][0][1]['x'][gl:-gr], self_conv_test_vect_pw[1][1],'.-')
    plt.grid()
    plt.ylim(ylim_inf,ylim_sup)
    plt.xlabel('space')
    #plt.ylabel('(h1-h2)-4(h2-h3)')
    data_name = "./data"+str(data[0])+"/"+field+"_epsilon"+str(epsilon[n_ind])+"dx"+str(dx)+"_pointwise_convergence.mp4"

    def animate1(i):
        line1.set_ydata(self_conv_test_vect_pw[i][0])  # update the data.
        line2.set_ydata(self_conv_test_vect_pw[i][1])
        return line1,line2

    #def animate2(i):
     #   line2.set_ydata(self_conv_test_vect_pw[i][1])  # update the data.
      #  return line2,


    ani = animation.FuncAnimation(
        fig, animate1, interval=80, blit=True, frames=number_steps-1)
    ani.save(data_name)

    
def animate_one_field(field_number,h_ind,big_DF_ind,model,data,big_DF,ylim_inf,ylim_sup,number_steps):
    fig, ax = plt.subplots()
    field = field_number
    data_name = "./data"+str(data[0])+"/field"+str(field)+".mp4"
    line, = ax.plot(big_DF[big_DF_ind][h_ind][0]['x'], big_DF[big_DF_ind][h_ind][0]['field'+str(field)])
    plt.ylim(ylim_inf,ylim_sup)
    plt.xlabel('space')
    plt.ylabel('field'+str(field))
    #plt.xlim(-5,0)
    plt.grid()

    def animate(i):
        #line.set_xdata(DF[0][i]['x'])
        line.set_ydata(big_DF[big_DF_ind][h_ind][i]['field'+str(field)])  # update the data.
        return line,


    ani = animation.FuncAnimation(
        fig, animate, interval=50, blit=True, save_count=number_steps-1)
    ani.save(data_name)
    
def animate_multiple_fields(big_DF,fields_to_print,h_ind,big_DF_ind,model,data,dx,gl,gr,ylim_inf,ylim_sup,number_steps,domain_lenght):
    data_name = "./data"+str(data[0])+"/fields.mp4"
    fig = plt.figure()
    ax1 = plt.axes(ylim=(ylim_inf,ylim_sup),xlim=(-float(dx)*gl,domain_lenght+float(dx)*gr))
    line, = ax1.plot([], [])
    plt.xlabel('x')
    plt.ylabel('fields')
    plotlays = [fields_to_print]
    lines = []
    for index in (fields_to_print):
        lobj = ax1.plot([],[],'-.',lw=3,alpha=0.6,label=index)[0]
        lines.append(lobj)

    x = []
    y = []
    for i in (fields_to_print):
        x.append([big_DF[0][h_ind][0]['x']]),
        y.append([ big_DF[0][h_ind][0][i]])


    def init():
        for line in lines:
            line.set_data(x,y)      
            #line.set_data([x1,x2,x3,x4,x5,x6,x7],[y1,y2,y3,y4,y5,y6,y7])
        return lines

    def animate(i):
        for n,j in enumerate (fields_to_print):
            x[n] = big_DF[big_DF_ind][h_ind][0]['x']
            y[n] = big_DF[big_DF_ind][h_ind][i][j]
        xlist = x
        ylist = y

        #for index in range(0,1):
        for lnum,line in enumerate(lines):
            line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately. 

        return lines

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=int(number_steps-1), interval=75, blit=True)
    plt.legend()
    plt.grid()
    plt.show()
    anim.save(data_name)
    
    
def theoretical_comparison(big_DF,theoretical_function,fields_to_print,h_ind,big_DF_ind,model,data,dx,gl,gr,number_steps,domain_lenght,ylim_inf,ylim_sup,times,amplitude_vector,ds,C):
    data_name = "./data"+str(data[0])+"/theoretical_comparison.mp4"
    fig = plt.figure()
    ax1 = plt.axes(ylim=(ylim_inf,ylim_sup),xlim=(-float(dx)*gl,domain_lenght+float(dx)*gr))
    line, = ax1.plot([], [])
    plt.xlabel('x')
    plt.ylabel('fields')
    plotlays = [fields_to_print]
    lines = []
    #t = time
    for index in (fields_to_print[0:-1]):
        lobj = ax1.plot([],[],'+',lw=3,alpha=0.6,label=index+"time"+str(t))[0]
        lines.append(lobj)
    lobj = ax1.plot([],[],lw=3,alpha=0.6,label="theoretical solution")[0]
    lines.append(lobj)

    x = []
    y = []
    for i in (fields_to_print[0:-1]):
        x.append([big_DF[big_DF_ind][h_ind][0]['x']]),
        y.append([ big_DF[big_DF_ind][h_ind][0][i]])

    x.append([big_DF[0][h_ind][0]['x']]),
    y.append([theoretical_function(big_DF[big_DF_ind][h_ind][0]['x'],times[0],float(amplitude_vector[big_DF_ind]),1,C,ds,domain_lenght)])
    #y.append([hyperbolic_chi_we_solution(big_DF[big_DF_ind][h_ind][0]['x'],times[0],1,0.5,domain_lenght)])

    def init():
        for line in lines:
            line.set_data(x,y)      
            #line.set_data([x1,x2,x3,x4,x5,x6,x7],[y1,y2,y3,y4,y5,y6,y7])
        return lines





    def animate(i):
        for n,j in enumerate (fields_to_print[0:-1]):
            x[n] = big_DF[big_DF_ind][h_ind][0]['x']
            y[n] = big_DF[big_DF_ind][h_ind][i][j]
        x[1] = big_DF[big_DF_ind][h_ind][0]['x']
        #y[1] = theoretical_function(big_DF[big_DF_ind][h_ind][0]['x'],times[i],float(amplitude_vector[big_DF_ind]))
        y[1] = theoretical_function(big_DF[big_DF_ind][h_ind][0]['x'],times[i],float(amplitude_vector[big_DF_ind]),1,C,ds,domain_lenght)
        t = times[i]

        xlist = x
        ylist = y

        #for index in range(0,1):
        for lnum,line in enumerate(lines):
            line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately. 

        return lines

    # call the animator.  blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=int(number_steps-1), interval=50, blit=True)
    plt.legend()
    plt.grid()
    plt.show()
    anim.save(data_name)
### scripts to call in the main ###
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
