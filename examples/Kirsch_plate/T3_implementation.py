# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:47:59 2020

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

Procedures = {"solver": {"type": "linear"}, 
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }

#---------------------------------------------------------------

# Material info
print("Defining material sets...")
MaterialSets = {

    '1': {'element': 'triangle',
          'plane deformation'    : 'plane stress',
          'material behavior': 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : 3*1e7, "poisson ratio" : 0.25},
          'geometric properties' : {'thikness': 0.5},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          
          }
    }

#---------------------------------------------------------------

# Mesh info
print("Meshing...")
mesh              = mesh_engine.Mesh()

mesh.Elements     = 2

mesh.NodesElement = 3
mesh.dofsNode     = 2
mesh.d            = 2

mesh.elements     = np.array([[0,1,2],
                              [0,2,3]])

mesh.elementMaterialTag = np.array([1,1])

mesh.elementType = np.array(["triangle",
                             "triangle"])
Lx = 3
Ly = 2
mesh.points  = np.array([[0.,  0., 0.],
                         [Lx, 0., 0.],
                         [Lx, Ly, 0.],
                         [0., Ly, 0.]])
mesh.Nodes        = len(mesh.points)

# Load info
print("Applying BCs...\n")

BCs = BC_engine.BoundaryConditions()

P = 1000

Dirichlet = 1
Neuman = 0
BCs.data = np.array([[0, Dirichlet, 0,  0],
                     [0, Dirichlet, 1,  0],
                     [1, Dirichlet, 1,  0],
                     [3, Dirichlet, 0,  0],
                     [3, Dirichlet, 1,  0],
                     [2, Neuman   , 1, -P]])#,
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

#post.run(mesh, U, MaterialSets, parameters)