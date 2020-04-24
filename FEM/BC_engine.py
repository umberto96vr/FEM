# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:50:27 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

import numpy as np

class BoundaryConditions():
    
    """ Class storing Von Neuman and Dirichlet boundary conditions
        
        Methods:
            self.add_load(node, load = None):
                Set the generalized force to the "node" numbered node to 0.
                
            self.constraint(node):
                Set the displacment to the "node" numbered node to 0.
    
    """
    

    
    def __init__(self):
        
        self.force  = None
        self.disp   = None
    
    def apply(self, K, F, dofs):
                
        if np.shape(self.force) == (2,):

            F[self.force[0]] += self.force[1]

        else:

            for (i, value) in enumerate(self.force[:,1]):

                loc = self.force[i,0]

                F[loc] = value
            

        if np.shape(self.disp) == (2,): 

            i   = self.disp[0]   

            K[:,loc]    = 0
            K[loc,:]    = 0
            K[loc,loc]  = 1            
            F[i]        = self.disp[1]

        else:

            for (i, value) in enumerate(self.disp[:,1]):

                loc = self.disp[i,0]

                K[:,loc]    = 0
                K[loc,:]    = 0
                K[loc,loc]  = 1
                F[loc]      = value
          
        
        return K, F
    
    def __str__(self):
        
        load_str = "Displaing boundary conditions...\n\n"
        load_str += "Constrained nodes:\n"
        for i in self.displacments:
            load_str += "{}, ".format(i)
        load_str += "\n\nLoad vector F: \n{}\n".format(self.F)    
        
        return load_str