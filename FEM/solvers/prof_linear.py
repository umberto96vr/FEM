#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:09:53 2020

@author: angelo.simone@unipd.it
"""

import numpy as np
from numpy.linalg import inv

import modules.FEM_engine as FEM_engine

# -----------------------------------------------------------------------------

def linear(mesh, BCs, MaterialSets, parameters):

    print("   . inizializing global arrays")
    K = np.zeros(shape=(mesh.dofs,mesh.dofs)) # allocate memory for global stiffness matrix
    F = np.zeros(shape=(mesh.dofs,1)) # allocate memory for global stiffness matrix

    print("   . generating global stiffness matrix")
    for e in range(mesh.Elements):
      
        k = FEM_engine.stiffness_matrix(e,mesh,MaterialSets,parameters)
        #print(k)
      
        # extract system dofs associated with element
        dof = FEM_engine.DofMap(e,mesh) 
        #print("Global to local mapping:",dof) 
    
        # assemble element matrices into global stiffness matrix
        K = FEM_engine.assemble(K,k,dof)

    print("   . applying boundary conditions")
    BCs.apply(K, F, mesh.dofs)

    print("   . solving KU=F for U")
    U = np.matmul(inv(K),F) # equivalent to inv(K).dot(F)

    return U
        