# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 19:59:59 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np
from numpy.linalg import inv
import time

import FEM.FEM_engine as FEM_engine

# -----------------------------------------------------------------------------

def linear(mesh, BCs, MaterialSets, parameters):

    K = np.zeros(shape = (mesh.dofs(),mesh.dofs()))
    F = np.zeros(shape = (mesh.dofs(),1))
    
    print("Assemblying global stiffness matrix...\n")
    
    start_a = time.time()
    
    for e in range(mesh.Elements):
    
        k = FEM_engine.stiffness_matrix(e, mesh, MaterialSets, parameters)
        
        # Get global dof associate with element e.
        dof = FEM_engine.DofMap(e,mesh)
        
        # Assemble the e-th local matrix into the global one.
        K = FEM_engine.assemble(K,k,dof)
    
    end_a = time.time()
    
    print("K:\n{}\n".format(K))
    print("Global stiffness matrix assembled in {}s\n".format(end_a - start_a))
    
    (Kr, F) = BCs.apply(K, F, mesh.dofs)
    
    print("Solving F = Ku...\n")
    U = np.matmul(inv(Kr),F)
    R = np.matmul(K, U)
    
    return U, R, K
        
    