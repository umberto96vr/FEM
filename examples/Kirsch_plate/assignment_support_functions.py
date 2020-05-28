# -*- coding: utf-8 -*-
"""
Created on Sun May 24 15:34:13 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np

"""
Infinite plate with circular hole uniformily loaded: Kirsch solution.
Below the functions giving the displacements and the stresses
for the solution are coded.
"""

def u_x(x, y, R, E, v, sig):
    
    E_ = E/(1-v**2)
    v_ = v/(1-v)
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    U_x = ( (1 + v_)/E_ )* sig * ( 
                (1/(1+v_))*r*np.cos(t) +
                (2/(1+v_))*(R**2 / r)*np.cos(t) +
                0.5*(R**2 / r)*np.cos(3*t) - 
                0.5 * (R**4 / r**3)*np.cos(3*t))
    if x == 0:
        U_x = 0
    
    return U_x

def u_y(x, y, R, E, v, sig):
    
    E_ = E/(1-v**2)
    v_ = v/(1-v)    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)

    U_y = ( (1 + v_)/E_ )* sig * ( 
                (-v_ /(1+v_))*r*np.sin(t) -
                ((1 - v_)/(1+v_))*(R**2 / r)*np.sin(t) +
                0.5*(R**2 / r)*np.sin(3*t) - 
                0.5 * (R**4 / r**3)*np.sin(3*t))
    
    if y == 0:
        U_y = 0
    
    return U_y

def s_x(x, y, R, sig): 
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_x = sig*(1 - (R**2/r**2)*(1.5*np.cos(2*t) + np.cos(4*t)) + 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_x

def s_y(x, y, R, sig):
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_y = sig*( -(R**2/r**2)*(0.5*np.cos(2*t) - np.cos(4*t)) - 
               1.5*(R**4/r**4)*np.cos(4*t))
    
    return S_y

def s_xy(x, y, R, sig):
    
    if x != 0:       
        t = np.arctan2(y,x) #(y/x)
    else:
        t = np.pi/2
    r = np.sqrt(x**2 + y**2)
    
    S_xy = sig*( -(R**2/r**2)*(0.5*np.sin(2*t) + np.sin(4*t)) + 
               1.5*(R**4/r**4)*np.sin(4*t))
    
    return S_xy

def change_mesh_size(filename, mesh_size):
    """
    This function is specifically designed for the a specific
    conventinal .geo file structure adopted in this assignment
    """ 
    with open(filename +'.geo', 'r') as file:
        # read a list of lines into data
        data = file.readlines()
        
    # changing mesh size
    data[1] = 'h = {};\n'.format(mesh_size)
    
    # and write everything back
    with open(filename +'.geo', 'w') as file:
        file.writelines( data )
