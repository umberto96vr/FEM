# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 14:03:21 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

import numpy as np

#-----------------------------------------------------------------------------

def mid_point_quad(n_interval, a, b, func):
    """
    Returns value of the finite integral of f(x)
    in dx from x = a to x = b.

    Parameters
    ----------
    n_interval : Number of descritization
                 intervals
    a :          Starting point of integration.
    b :          Ending point of integration.
    func :       function func(x)

    Returns
    -------
    I:           Value of the integral

    """
    
    delta = (b-a)/n_interval
    I = 0
    xi = a + delta/2
    
    for i in range (n_interval):
         
        I += func(xi)*delta
        xi += delta
        
    return I  

#-----------------------------------------------------------------------------

# Find the number of integration points such as the absolute
# errore err is below 0.01  
# err = 1
# i = 0
# while err > 0.01:
#     print("i = {}".format(i))
#     i += 1 
#     I = mid_point_quad(i, 0, 3, func)
    
#     err = abs(I - 12)
#     print(I)
#     print("err(i) = {}\n".format(err))

#-----------------------------------------------------------------------------


def gauss_leg_quadr(n, a, b, func):
    """
    Gauss-Legendre numerical integration of f(x) in dx from a to b

    Parameters
    ----------
    n : Number of integration poits.
    
    a,b: domain boundaries along x. 
    
    func : f(x)

    Returns
    -------
    I : Integral value.

    """
    
    if n > 4:
        print("Error: maximum number of integration points is n = 4!")
        
        return 0
    
    else:
        
        # Quadrature scheme
        # quad_scheme_GL = {
        #                 2: np.array([
        #                             [-1/np.sqrt(3),1/np.sqrt(3)],
        #                             [1,1]]),
        #                 3: np.array([
        #                              [-np.sqrt(3/5), 0, np.sqrt(3/5)],
        #                              [5/9, 8/9, 5/9]]), 
        #                 4: np.array([
        #                             [-np.sqrt((15 + 2*np.sqrt(30))/35), 
        #                              -np.sqrt((15 - 2*np.sqrt(30))/35), 
        #                               np.sqrt((15 - 2*np.sqrt(30))/35), 
        #                               np.sqrt((15 + 2*np.sqrt(30))/35)],
        #                              [(18 - np.sqrt(30))/36,
        #                               (18 + np.sqrt(30))/36,
        #                               (18 + np.sqrt(30))/36,
        #                               (18 - np.sqrt(30))/36]])}
        quad_scheme_GL = {
                        2: np.array([
                                    [-0.5773502691896258,
                                      0.5773502691896258],
                                    [1,1]]),
                        3: np.array([
                                     [-0.7745966692414834, 
                                       0, 
                                       0.7745966692414834],
                                     [0.5555555555555556, 
                                      0.8888888888888888, 
                                      0.5555555555555556]]), 
                        4: np.array([
                                    [-0.8611363115940526, 
                                     -0.33998104358485626, 
                                      0.33998104358485626, 
                                      0.8611363115940526],
                                     [0.34785484513745385,
                                      0.6521451548625462,
                                      0.6521451548625462,
                                      0.34785484513745385]])}
        
        I = 0   # Initialize integral value
        
        # Computation loop
        
        for i in range(n):
            
            xi = quad_scheme_GL[n][0][i]
            wi = quad_scheme_GL[n][1][i]
            c = (b-a)/2
            d = (b+a)/2
            
            I +=  wi*func(xi*c + d)
        
        I *= (b-a)/2
        
        return I

#-----------------------------------------------------------------------------

def gauss_leg_quadr2d(n, a, b, c, d, func):
     
    """
    Gauss-Legendre numerical integration of f(x,y) in dx and dy 
    for a rectangular domain spanning from a to b along x and from 
    c to d along y

    Parameters
    ----------
    n:    Number of integration poits.
    
    a,b:  domain limits along x.
    
    c,d:  domain limits along y.
    
    func: f(x,y)

    Returns
    -------
    I:    Integral value.

    """
    
    if n > 4:
        # n > 4 currently WIP
        print("Error: maximum number of integration points is n = 4!")
        return 0
    else:
        
        # Change of integration variable constants:
            # x = c1*xi  + c2
            # y = c3*eta + c4
            
        c1 = (b-a)/2
        c2 = (b+a)/2
        c3 = (d-c)/2
        c4 = (d+c)/2
            
        #Jacobian
        J = c1*c3   
    
        # Quadrature scheme
        quad_scheme_GL = {
                        2: np.array([
                                    [-0.5773502691896258,
                                      0.5773502691896258],
                                    [1,1]]),
                        3: np.array([
                                     [-0.7745966692414834, 
                                       0, 
                                       0.7745966692414834],
                                     [0.5555555555555556, 
                                      0.8888888888888888, 
                                      0.5555555555555556]]), 
                        4: np.array([
                                    [-0.8611363115940526, 
                                     -0.33998104358485626, 
                                      0.33998104358485626, 
                                      0.8611363115940526],
                                     [0.34785484513745385,
                                      0.6521451548625462,
                                      0.6521451548625462,
                                      0.34785484513745385]])}
        
        I = 0   # Initialize integral value
        
        # Computation loop
        
        for i in range(n):
            
            xi_i = quad_scheme_GL[n][0][i]
            w_i = quad_scheme_GL[n][1][i]
            
            x_i = c1*xi_i + c2
            
            
            for j in range(n):
            
                eta_j = quad_scheme_GL[n][0][j]
                w_j = quad_scheme_GL[n][1][j]
            
                y_j = c3*eta_j + c4
                    
                I +=  w_i*w_j*func(x_i, y_j)
        
        I *= J
        
        return I


def gauss_leg_quadr2d_tri(n, X, Y, func):
     
    """
    Gauss-Legendre numerical integration of f(x,y) in dx and dy 
    for triangular domain of coordinates specified by vectors X and Y

    Parameters
    ----------
    n:    Number of integration poits.
    
    X:    X = np.array([x1, x2, x3])
    
    Y:    Y = np.array([y1, y2, y3])
    
    func: f(x,y)

    Returns
    -------
    I:    Integral value.

    """
    
    
    if n == 1 or n == 3:
        
        
        # Change of integratio variable constants:
        def xx(xi,eta):
            
            return X[0] + (X[1]-X[0])*xi + (X[2]-X[0])*eta
        
        def yy(xi,eta):
            
            return Y[0] + (Y[1]-Y[0])*xi + (Y[2]-Y[0])*eta
            
        #Jacobian
        J = (X[1]-X[0])*(Y[2]-Y[0]) - (X[2]-X[0])*(Y[1]-Y[0])
    
        # Quadrature scheme
        quad_scheme_GL = {
                        1: np.array([
                                    [1],
                                    [0.3333333333333333],
                                    [0.3333333333333333]]),
                        3: np.array([
                                     [0.3333333333333333,  
                                      0.3333333333333333,  
                                      0.3333333333333333],
                                     [0, 0.5, 0.5],
                                     [0.5, 0, 0.5]])}
        
        I = 0   # Initialize integral value
        
        # Computation loop
        
        for i in range(n):
            
            xi_i = quad_scheme_GL[n][1][i]
            eta_i = quad_scheme_GL[n][2][i]
            w_i = quad_scheme_GL[n][0][i]
            
            I +=  0.5*w_i*func(xx(xi_i, eta_i),yy(xi_i, eta_i))
        
        I *= J
        
        return I
    else:
        # n > 4 currently WIP
        print("Error: n has to be either 1 or 3!")
        return 0