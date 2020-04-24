# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 17:43:32 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import numpy as np


class QuadratureRuleError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class InversionError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class FemEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

class Integration():
    
    quad_scheme = {
        
        'Gauss-Legendre': { 'line' : {
                                    1: np.array([
                                                [2.],
                                                [0.]]),
                                    2: np.array([
                                                [-0.5773502691896258,
                                                 0.5773502691896258],
                                                [1.,1.]]),
                                    3: np.array([
                                                [-0.7745966692414834, 
                                                 0., 
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
                                                 0.34785484513745385]])},
                    'rectangular' : {
                                    1: np.array([
                                                [2.],
                                                [0.]]),
                                    2: np.array([
                                                [-0.5773502691896258,
                                                 0.5773502691896258],
                                                [1.,1.]]),
                                    3: np.array([
                                                [-0.7745966692414834, 
                                                 0., 
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
                                                 0.34785484513745385]])},
                     'triangular' : {
                                    1: np.array([
                                                [1.],
                                                [0.3333333333333333],
                                                [0.3333333333333333]]),
                                    3: np.array([
                                                [0.3333333333333333,  
                                                 0.3333333333333333,  
                                                 0.3333333333333333],
                                                [0., 0.5, 0.5],
                                                [0.5, 0., 0.5]])}
                             }
        }
        
    
    
    def __init__(self):
        pass
    
    def quadrature_rule(self, rule, domain, n_points):
        
        if rule in self.quad_scheme.keys():
            
            if domain in self.quad_scheme[rule].keys():
                
                if n_points in self.quad_scheme[rule][domain].keys():
                    
                    weights = self.quad_scheme[rule][domain][n_points][1,:]
                    points  = self.quad_scheme[rule][domain][n_points][0,:]
                    
                    return (weights, points)
                else:
                    
                    err_msg = "'{}' quadrature for a ".format(rule)
                    err_msg +="'{}' domain hasn't been defined for ".format(domain)
                    err_msg+= "{} points of integration!".format(n_points)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()
                    
            else:
                    
                    err_msg = "'{}' quadrature for a ".format(rule)
                    err_msg +="'{}' domain hasn't been defined!".format(domain)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()
        else:
                    
                    err_msg = "'{}' quadrature hasn't been defined!".format(rule)
                    
                    raise QuadratureRuleError(err_msg)
                    
                    # print('Error: ', err_msg)
                    # sys.exit()

def shape_functions(elem_coordinates, eval_point, NodesElement, d):
    """
    Return the shape function derivatives vector evalueted
    at the quadrature scheme integration point. \n
    The conversion from the [-1, 1] domain to the [x_i, x_j] domain
    is done internally in this function

    Parameters
    ----------
    ElemCoordinates : Element coordinates xi, xj for a certain element \n
    eval_point : Shape function derivatives vector evaluation point in [-1, 1] domain.\n
    NodesElement: number of nodes in the element.\n
    d: spatial dimension of the problem, e.g. 1 for 1D, 2 for 2D,...\n

    Returns
    -------
    dN: shape function vector

    """
    
    x = elem_coordinates.reshape(NodesElement,d)
    
    # Vector initialisation
    N  = np.zeros( [d, NodesElement] )
    dN = np.zeros( [d, NodesElement] )
    J  = np.zeros( [d,d] )
    
    xi = eval_point
    
    # Shape functions evaluation.
    N[0,0] = -0.5*xi + 0.5
    N[0,1] = 0.5*xi + 0.5
    
    # Shape functions derivatives evaluation
    dN[0,0] = -0.5
    dN[0,1] =  0.5
    
    # Jacobian matrix evaluation
    J = np.dot(dN,x)
    
    if d <= 2:
        
        detJ = determinant(J)
        invJ = inverse(J)
    else:
        invJ = np.linalg.inv(J)
        detJ = np.linalg.det(J)
    
    if detJ < 0:
        raise FemEngineError("The Jacobian |J| < 0!")
    else:    
        
        # dN has to be expressed in terms of the global coordinate system (x,y,z)
        # The integrand function is F(x,y)*|J|*d_xi * d_eta that is
        #   B(x,y)^T * D * B(x,y) * |J| * d_xi * d_eta 
        # = B(xi,eta)^T * (J^-1)^T* D * B(xi,eta) * J^T * |J| * d_xi * d_eta 
            
        dN = np.dot(invJ, dN)
        
        return (detJ, dN)               


        
def determinant(matrix):
    """
    Compute the determinant of the given matrix.
    The matrix 

    Parameters
    ----------
    matrix 

    Returns
    -------
    det: Determinant

    """   
    det = 0
    n = np.shape(matrix)[1]
    m = np.shape(matrix)[1]
    
    if n == m:        
    
        if n == 1:
            return matrix[0,0]
        else:
            
            for i in range(n):
                minor = np.array( [ matrix[1:n, j] for j in range(m) if j != i] )
                det += (-1)**(i) * matrix[0,i] * determinant(minor)
            
            return det
            
    else:
        raise InversionError("The matrix isn't square!")
                

def inverse(matrix):
    """
    Compute the determinant of the given matrix.
    The matrix 

    Parameters
    ----------
    matrix 

    Returns
    -------
    det: Determinant

    """   
    det = 0
    n = np.shape(matrix)[1]
    m = np.shape(matrix)[1]
    
    if n == m:        
    
        if n == 1:
            return 1/matrix[0,0]
        elif n == 2:
            det = determinant(matrix)
            
            return (1/det) * np.array( [[ matrix[1,1], -matrix[0,1] ],
                                        [-matrix[1,0],  matrix[0,0] ]] )
            
            
        else:
            raise InversionError("Code for matrix bigger than 3x3 hasn't been developed yet!")
            
    else:
        raise InversionError("The matrix isn't square!")         
    
        
        