# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:50:27 2020

@author: kaise
"""

import numpy as np

class BoundaryConditions():
    
    """ Class storing Von Neuman and Dirichlet boundary conditions"""

    
    def __init__(self, loads = None, constraints = None):
        
        self.loads = loads
        self.constraints = constraints
    
    