# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:38:10 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import time
sys.path.append("../../source")

import numpy as np

from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
import FEM.post.postprocess as post

np.set_printoptions(precision=4)

#---------------------------------------------------------------
# Clearing variable workspace

#from IPython import get_ipython;   
#get_ipython().magic('reset -sf')

#---------------------------------------------------------------
# Parameters
print("\n--> Pre-processing...\n")
print("Initializing model parameters...")

Procedures = {"solver": {"type": "linear"}, # try with nonlinear
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }

#---------------------------------------------------------------

# Material info
print("Defining material sets...")
MaterialSets = {

    '1': {'element': 'bar2d',
          'elastic properties'   : {"Young modulus" : 206000, "poisson ratio" : 0.3},
          'geometric properties' : {'area': 500},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          }
    }

#---------------------------------------------------------------

# Mesh info
print("Meshing...")
mesh              = mesh_engine.Mesh()

mesh.Elements     = 15
mesh.NodesElement = 2
mesh.dofsNode     = 2
mesh.d            = 2

mesh.elements     = np.array([[0,5],
                              [0,1],
                              [5,1],
                              [5,6],
                              [1,6],
                              [1,2],
                              [6,2],
                              [6,7],
                              [2,7],
                              [2,3],
                              [7,3],
                              [7,8],
                              [3,8],
                              [3,4],
                              [8,4]])

mesh.elementMaterialTag = np.ones(len(mesh.elements), dtype = int)

mesh.elementType = np.array([["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"],
                             ["bar2d"]])

L1 = 4000
L2 = 6000
mesh.points  = np.array([[     0,  0],
                         [    L1,  0],
                         [  2*L1,  0],
                         [  3*L1,  0],
                         [  4*L1,  0],
                         [0.5*L1, L2],
                         [1.5*L1, L2],
                         [2.5*L1, L2],
                         [3.5*L1, L2]])


mesh.Nodes        = len(mesh.points)

# Load info
print("Applying BCs...\n")

BCs = BC_engine.BoundaryConditions()

q1 = -10*L1

Dirichlet = 1
Neuman = 0
BCs.data = np.array([[0, Dirichlet, 0,  0],
                     [0, Dirichlet, 1,  0],
                     [4, Dirichlet, 1,  0],
                     [1, Neuman   , 1, q1],
                     [2, Neuman   , 1, q1],
                     [3, Neuman   , 1, q1]])#,
                     #[4, Neuman   , 1, q1],
                     #[0, Neuman   , 1, q1]])

#---------------------------------------------------------------

print("--> Solving...\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)

end = time.time()


print("Total time: {}s".format(end-start))

print("\n-- Post-processing...\n")

U_print = U.copy()

print("   . solution U:\n{}\n".format(U_print.reshape(mesh.Nodes,2)))
print("   . reaction forces R:\n{}\n".format(R.reshape(mesh.Nodes,2)))

post.run(mesh, U, MaterialSets, Procedures)