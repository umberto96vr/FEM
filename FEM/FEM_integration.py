# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 17:43:32 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import numpy as np
import sys

class QuadratureRuleError(Exception):
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
                    
                    weights = self.quad_scheme[rule][domain][n_points][0,:]
                    points  = self.quad_scheme[rule][domain][n_points][1,:]
                    
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
                    
                    
        
        
        