# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 15:51:53 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

# Clearing variable workspace
#from IPython import get_ipython;   
#get_ipython().magic('reset -sf')

#---------------------------------------------------------------
import sys
import time
sys.path.append("../source")

import numpy as np

from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver

#---------------------------------------------------------------

# Parameters
print("Initializing model parameters...\n")

parameters = {"solver": {"type": "linear"}}
L = 3
P = 0.3
a = 1
n = 40
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
          'elastic properties'   : {"Young modulus" : 1, "poisson ratio" : 0},
          'geometric properties' : {'t': 1},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          },
    '7': {'element': 'quad',
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

mesh.Elements     = n
mesh.Nodes        = int(2+ 2*mesh.Elements)
mesh.NodesElement = 4
mesh.dofsNode     = 2
mesh.d            = 2


conect = np.zeros([mesh.Elements, 2+mesh.NodesElement], dtype = int)

for i in range(mesh.Elements):
    
    conect[i,:] = [4, 6, i, (i+1), (mesh.Elements+i+2), (mesh.Elements+i+1)]
    
mesh.elements = conect

coord = np.zeros([mesh.Nodes, mesh.d])

for i in range(mesh.Nodes):
    if i < int(mesh.Nodes/2):
        if i == 0:        
            coord[i,0] = 0        
        else:
            coord[i,0] = coord[i-1,0] + L/mesh.Elements
        
        coord[i,1] = 0
    else:
        coord[i,0] = coord[int(i - mesh.Nodes/2),0]
        coord[i,1] = a
mesh.coordinates = coord

# Load info
print("Applying BCs...\n")

BCs       = BC_engine.BoundaryConditions()
BCs.data = np.array([[0, 1, 0, 0],
                     [0, 1, 1, 0],
                     [mesh.Nodes-1, 1, 0, 0],
                     [mesh.Nodes-1, 0, 0, P/2],
                     [n, 0, 0, P/2]])

#---------------------------------------------------------------

print("-------------------------------------------------------")
print("                     SOLVING")
print("-------------------------------------------------------\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, parameters)

end = time.time()

for i,Ui in enumerate(U):
  if Ui < 1e-10:
    U[i] = 0

print("U: \n{}\n".format(U.reshape((len(U),1))))
print("R: \n{}\n".format(R.reshape((len(R),1))))

print("Total time: {}s".format(end-start))

X_old = mesh.coordinates
X_new = np.zeros(np.shape(X_old))
U_res = X_new.copy()

for k in range(len(U)):
    i = int((k - k%2)/2)
    j = k%2
    X_new[i,j] = X_old[i,j] + U[k]
    U_res[i,j] = U[k]
    
import matplotlib.pyplot as plt

plt.figure(1)
plt.plot(X_old[:,0],X_old[:,1], '*b')
plt.plot(X_new[:,0],X_new[:,1], '*r')


plt.figure(3)
plt.spy(K)

plt.show()

new_cord = np.zeros((mesh.Elements, mesh.NodesElement * mesh.d))
old_cord = new_cord.copy()

for i in range(mesh.Elements):
    for j,n in enumerate(mesh.nodesInElement(i)):
        old_cord[i, (2*j):(2*j+2)] = mesh.elem_coordinates(i)[j]
        new_cord[i, (2*j):(2*j+2)] =  old_cord[i, (2*j):(2*j+2)] + np.array([float(U[2*n]), float(U[2*n+1])])

