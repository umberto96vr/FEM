# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:38:10 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import time
sys.path.append("../source")

import numpy as np

from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver

#---------------------------------------------------------------
# Clearing variable workspace

#from IPython import get_ipython;   
#get_ipython().magic('reset -sf')

#---------------------------------------------------------------
# Parameters
print("Initializing model parameters...\n")

parameters = {"solver": {"type": "linear"}}

#---------------------------------------------------------------

# Material info
print("Defining material sets...\n")
MaterialSets = {

    '1': {'element': 'spring',
          'elastic properties'   : {'stiffness' : 2},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          },
    '2': {'element': 'bar',
          'elastic properties'   : {"Young modulus" : 2, "poisson ratio" : 0.3},
          'geometric properties' : {'area': 2},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          },
    '3': {'element': 'bar',
          'elastic properties'   : {"Young modulus" : 2, "poisson ratio" : 0.3},
          'geometric properties' : {'area': 1},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          },
    '4': {'element': 'bar',
          'elastic properties'   : {"Young modulus" : 2, "poisson ratio" : 0.3},
          'geometric properties' : {'area': 1},
          'stiffness matrix'     : {'evaluation': 'numerical integration',
                                    'domain'    : 'line',
                                    'rule'      : 'Gauss-Legendre',
                                    'points'    : 4}
          },
    '5': {'element': 'bar2d',
          'elastic properties'   : {"Young modulus" : 206000, "poisson ratio" : 0.3},
          'geometric properties' : {'area': 220},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          },
    '6': {'element': 'quad',
          'elastic properties'   : {"Young modulus" : 206000, "poisson ratio" : 0.3},
          'geometric properties' : {'t': 220},
          'stiffness matrix'     : {'evaluation': 'numerical integration',
                                    'domain'    : '2d',
                                    'rule'      : 'Gauss-Legendre',
                                    'points'    : 4}
          }
    }

#---------------------------------------------------------------

# Mesh info
print("Meshing...\n")
mesh              = mesh_engine.Mesh()

mesh.Elements     = 3
mesh.Nodes        = 4
mesh.NodesElement = 2
mesh.dofsNode     = 2
mesh.d            = 2

mesh.elements     = np.array([[3,5,0,1],
                              [3,5,1,3],
                              [3,5,1,2]])
L = 2000
mesh.coordinates  = np.array([[0  , 0  ],
                              [L  , L  ],
                              [2*L, L  ],
                              [L  , 2*L]])

# Load info
print("Applying BCs...\n")

BCs       = BC_engine.BoundaryConditions()
BCs.force = np.array([[3, -21560]])
BCs.disp  = np.array([[0,0],
                      [1,0],
                      [4,0],
                      [5,0],
                      [6,0],
                      [7,0]])
BCs.data = np.array([[]])

#---------------------------------------------------------------

print("-------------------------------------------------------")
print("                     SOLVING")
print("-------------------------------------------------------\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, parameters)

end = time.time()

print("U: \n{}\n".format(U.reshape((len(U),1))))
print("R: \n{}\n".format(R.reshape((len(R),1))))

print("Total time: {}s".format(end-start))