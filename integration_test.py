# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 17:18:40 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import numpy as np
from FEM import integration as inte

#-----------------------------------------------------------------------------
#                       TESTING FUNCTIONS
#-----------------------------------------------------------------------------

def func(x):
    
    return x**2 +1

#-----------------------------------------------------------------------------

def func2d(x,y):
    
    return (x**2)*(y**2)

#-----------------------------------------------------------------------------

def circ(x):
    return np.sqrt(1 - x**2)

#-----------------------------------------------------------------------------

print("\n#-----------------------------------------------------------------------------")
print("#                               TESTING")
print("#-----------------------------------------------------------------------------\n")
print("Integral I of f(x) = x^2 +1 for x from 0 to 3\n\n")

print("Mid-point quadrature, n_interval = 13")

print("I = {}\n".format(inte.mid_point_quad(13, 0, 3, func)))

print("Gauss-Legendre quadrature, n = 4")

print("I = {}\n".format(inte.gauss_leg_quadr(4, 0, 3, func))) 

print("\n-----------------------------------------------------------------------------\n")

print("Integral I of f(x,y) = x^5 * y^5 for x,y both from 0 to 2\nusing Gauss-Legender with n = 4 integration points\n")

print("I = {}\n".format(inte.gauss_leg_quadr2d(4, 0, 2, 0, 2, func2d)))
    
print("\n#-----------------------------------------------------------------------------")
print("#                               COMPUTING PI")
print("#-----------------------------------------------------------------------------\n")   
    
    
for j in [2,3,4]:
    
    I = 4*inte.gauss_leg_quadr(j, 0, 1, circ)
    print('Gauss-Legendre, n = {}'.format(j))
    print('pi = {}'.format(I))
    print("err% = {}%\n".format((I-np.pi)/(0.01*np.pi)))
    
for j in [10, 30, 50, 70, 90]:
    
    I = 4*inte.mid_point_quad(j, 0, 1, circ)
    print('Mid-point quadrature, n = {}'.format(j))
    print('pi = {}'.format(I))
    print("err% = {}%\n".format((I-np.pi)/(0.01*np.pi)))
    
print("\n#-----------------------------------------------------------------------------")
print("#                               TESTING TRIANGULAR")
print("#-----------------------------------------------------------------------------\n")

X = np.array([0, 1, 1])
Y = np.array([0, 1, 2])
def test_tri(x,y):
    
    return x**2 + y**2

print("Integral I of f(x) = sqrt(x + y) over the unitary triangular domain\n\n")

print("Gauss-Legendre quadrature, n = 1")

print("I = {}\n".format(inte.gauss_leg_quadr2d_tri(1, X, Y, test_tri)))

print("Gauss-Legendre quadrature, n = 3")

print("I = {}\n".format(inte.gauss_leg_quadr2d_tri(3, X, Y, test_tri))) 

