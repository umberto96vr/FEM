# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:50:27 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

#import numpy as np

class BoundaryConditions():
    
    """ Class storing Von Neuman and Dirichlet boundary conditions
        
        Methods:
            self.add_load(node, load = None):
                Set the generalized force to the "node" numbered node to 0.
                
            self.constraint(node):
                Set the displacment to the "node" numbered node to 0.
    
    """
    

    
    def __init__(self, loads = None, constraints = None):
        
        self.loads = loads
        self.constraints = constraints
    
    def add_load(self, node, load = None):
        
        self.loads[node] = load
    
    def __str__(self):
        
        load_str = "Displaing boundary conditions...\n\n"
        load_str += "Constrained nodes:\n"
        for i in self.constraints:
            load_str += "{}, ".format(i)
        load_str += "\n\nLoad vector F: \n{}\n".format(self.loads)    
        
        return load_str