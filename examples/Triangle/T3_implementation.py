# -*- coding: utf-8 -*-
"""
Created on Thu May 14 16:43:28 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

import sys
import time
sys.path.append("../../source")

import numpy as np
import meshio
from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
import FEM.post.post2d as post2d

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

    'box': {'element': 'triangle',
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

mesh_file = "square1"
mesh = mesh_engine.GMSH(mesh_file)


# Load info
print("Applying BCs...\n")

BCs = BC_engine.BoundaryConditions()

P = 1000.

Dirichlet = 1
Neumann = 0
BC1 = BCs.set_bc("Dirichlet", "Fixed", mesh, [0, 1], [.0, .0])
BC2 = BCs.set_bc("Neumann", "Applied vertical force", mesh, [1], [-P])
BC3 = BCs.set_bc("Dirichlet", "Horizontal roller", mesh, [1], [.0])

BCs.data = BC1 + BC2 + BC3

#---------------------------------------------------------------

print("--> Solving...\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)

end = time.time()


print("Total time: {}s".format(end-start))

print("\n-- Post-processing...\n")

U_print = U.copy().reshape(len(mesh.points),mesh.d)

print("   . solution U:\n{}\n".format(U_print))
print("   . reaction forces R:\n{}\n".format(R.reshape(len(mesh.points),mesh.d)))

(stress, von_mises) = post2d.run_post(mesh, MaterialSets, U_print)
