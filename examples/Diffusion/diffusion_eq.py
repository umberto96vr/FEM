# -*- coding: utf-8 -*-
"""
Created on Tue May  5 11:31:36 2020

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
from FEM.FEM_utilities import element_table
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#---------------------------------------------------------------
# Clearing variable workspace

#from IPython import get_ipython;   
#get_ipython().magic('reset -sf')

#---------------------------------------------------------------
# Parameters
print("Initializing model parameters...\n")

dt = (60*60*24*365)*1000
n_time = 5000

time_tot = n_time*dt

Procedures = {"solver": {"type":'nonlinear', 
                          "dt": dt,
                          "Tot time": time_tot}, # try with nonlinear
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }



#---------------------------------------------------------------

# Material info
print("Defining material sets...\n")
MaterialSets = {

    '1': {'element': 'diff_bar',
          'elastic properties'   : {'diffusivity' : (11/9)*1e-06},
          'mass properties'      : {'density'     : 1},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          }
    }

#---------------------------------------------------------------

# Mesh info
print("Meshing...\n")

L = 10000
n = 40

mesh              = mesh_engine.Mesh()

mesh.Elements     = n
mesh.Nodes        = n+1
mesh.NodesElement = 2
mesh.dofsNode     = 1
mesh.d            = 1
mesh.points       = np.zeros(shape = (mesh.Nodes,3))
mesh.points[:,0]  = np.linspace(0,L,mesh.Nodes)
mesh.elements     = element_table(mesh)
mesh.elementMaterialTag = np.ones(n, dtype = int)
mesh.elementType = np.array([["diff_bar"] for i in range(n)])

# Load info
print("Applying BCs...\n")

BCs       = BC_engine.BoundaryConditions()

F_size = mesh.Nodes

data = np.zeros([2, 4])

S = 9.6296*1e-13

data[0,:] = np.array([0,1,0,0]) 
data[1,:] = np.array([-1,1,0,0]) 

BCs.data = data
BCs.init_cond = np.zeros(F_size)
#---------------------------------------------------------------

print("-------------------------------------------------------")
print("                     SOLVING")
print("-------------------------------------------------------\n")

# Global solver time assesment
start = time.time()

U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)

end = time.time()

#print("U: \n{}\n".format(U.reshape((len(U),1))))
#print("R: \n{}\n".format(R.reshape((len(R),1))))

print("Total time: {}s".format(end-start))

# Result plot

tempi = [50, 100, 150]
plt.figure(1)
plt.grid(which = 'both')
for t in tempi:
    plt.plot(mesh.points[:,0], U[:,t], '*', markersize = 3)
plt.ylabel(r'Temperature (°C)')
plt.xlabel(r'Distance (m)')
plt.legend(['t = {}000 years'.format(i) for i in tempi])

