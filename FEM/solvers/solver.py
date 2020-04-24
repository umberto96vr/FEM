# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 19:51:19 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)


def run(mesh, BCs, material, parameters):
    
    possible_solvers = ['linear', 'nonlinear']

    solverType = parameters["solver"]["type"]
    
    if solverType in possible_solvers:
    
        exec("from FEM.solvers."+solverType+" import "+solverType)
        
        (U, R, K) = eval(solverType+"(mesh, BCs, material, parameters)")    
        
        # ...you might want to have a safer approach to eval:
        # https://www.geeksforgeeks.org/eval-in-python/
        
        return (U, R, K)
    else:
        raise SolverError("Invalid solver type!")

