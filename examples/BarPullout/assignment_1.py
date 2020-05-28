# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 19:37:38 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

#---------------------------------------------------------------
import sys
import time
sys.path.append("../../source")

import numpy as np

from FEM import mesh_engine 
from FEM import BC_engine
from FEM.FEM_utilities import element_table
import FEM.solvers.solver as solver
import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('figure', max_open_warning = 1000)

elem_plot = True
conv_plot = True
save_elem = False
save_conv = False

#---------------------------------------------------------------
# Output path and parameters
out_path_eps = "../../../Report/Figures_eps/"
out_path_pdf = "../../../Report/Figures_pdf/"
out_path_png = "../../../Report/Figures_png/"

markers = ['*','^','x','d','o']

#---------------------------------------------------------------
# Parameters
print("Initializing model parameters...\n")

parameters = {"strain components": 1, #TODO: to be related to the spatial dimension
              #strain component -> stress/strain state: 1D, 2D plane stress, 2D plane strain, 3D
              
              "solver": {"type": "linear"}, # try with nonlinear
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }
L = 10
P = 1
E = 1
A = 1
alfa_range = [0.5, 1, 2, 4]
k_range = [E*A*i**2 for i in alfa_range]
n_cases = [2**i for i in range(2,7)]
h_sizes = [L/n for n in n_cases]

ERR = np.zeros((len(alfa_range), len(n_cases)))
#---------------------------------------------------------------

# Plottin indexex parameters
k_max = 2
i_max = len(alfa_range)
j_max = len(n_cases)

#---------------------------------------------------------------

# Error analysis
e = np.zeros(shape = (k_max, i_max, j_max))
eps_tab = e.copy()
tab = e.copy()
eps = [((alfa_range[i] * P**2)/(2*k_range[i]))*(1 - np.exp(-2*alfa_range[i]*L)) for i in range(len(k_range))]

#---------------------------------------------------------------


for (k, NodesElement) in enumerate([2,3]):
    for (i, alfa) in enumerate(alfa_range):
        # Material info
        print("Defining material sets...\n")
        MaterialSets = {       
            
            '1': {'element': 'bar_pullout',
                  'elastic properties'   : {"Young modulus" : E, 
                                            "poisson ratio" : 0.3,
                                            "k" : k_range[i]},
                  'geometric properties' : {'area': A},
                  'stiffness matrix'     : {'evaluation': 'numerical integration',
                                            'domain'    : 'line',
                                            'rule'      : 'Gauss-Legendre',
                                            'points'    : 3}
                  },
            '2': {'element': 'bar_pullout',
                  'elastic properties'   : {"Young modulus" : E, 
                                            "poisson ratio" : 0.3,
                                            "k" : k_range[i]},
                  'geometric properties' : {'area': A},
                  'stiffness matrix'     : {'evaluation': 'closed form'}
                  }
            }
        
        #---------------------------------------------------------------
        for (j,n) in enumerate(n_cases):
            # Mesh info
            print("Meshing...\n")
            mesh              = mesh_engine.Mesh()
            
            mesh.Elements     = n            
            mesh.NodesElement = NodesElement
            mesh.points       = np.zeros(shape = (int((mesh.NodesElement-1)*mesh.Elements+1),3))
            mesh.points[:,0]  = np.linspace(0, L, len(mesh.points[:,0]))
            mesh.Nodes        = (mesh.NodesElement-1)*mesh.Elements+1
            mesh.dofsNode     = 1
            mesh.d            = 1
            #mesh.dofs         = mesh.Nodes*mesh.dofsNode
            mesh.elements     = element_table(mesh)
            mesh.elementMaterialTag = np.ones(n, dtype = int)
            mesh.elementType = np.array([["bar_pullout"] for i in range(n)])
            
            
            # Load info
            print("Applying BCs...\n")
            
            BCs       = BC_engine.BoundaryConditions()
            BCs.data  = np.array([[0, 0, 0, -P],
                                  [int(mesh.Nodes-1),
                                   1, 0,
                                   -P/(E*A*alfa) * np.exp(-alfa*L)]])
    
            
            print("-------------------------------------------------------")
            print("                     SOLVING")
            print("-------------------------------------------------------\n")
            
            # Solving
            (U, R, K) = solver.run(mesh, BCs, MaterialSets, parameters)
            
            print("------------------------ END --------------------------\n")
            
            # Error analysis
            eps_h = float(0.5 * U.T @ K @ U  )
            
            e[k,i,j] = np.sqrt((eps[i] - eps_h) / eps[i])
            
            # Plotting
            if elem_plot == True:
                plt.figure(j+ k*i_max*j_max + i*j_max + 1)
                plt.plot(mesh.points[:,0], U, '-xb')
                true_sol = -P/(E*A*alfa) * np.exp(-alfa*mesh.points[:,0])
                plt.plot(mesh.points[:,0], true_sol, '-+r')
                plt.title("Solution comparison: {}-node bar element,\n{} elements".format(NodesElement, n))
                plt.xlabel(r"x")
                plt.ylabel(r"u(x)")
                plt.legend(["Numerical solution", "Analytical solution"])
                plt.grid()
                if mesh.elements[0,1] == 1:
                    filename = "alfa{}_{}elems_{}-node_isoparametric".format(alfa,n,NodesElement)
                elif mesh.elements[0,1] == 2:
                    filename = "alfa{}_{}elems_{}-node_closedform".format(alfa,n,NodesElement)
                if save_elem == True:
                    plt.savefig(out_path_eps + "Elements/alfa{}/".format(alfa) + filename + ".eps")
                    plt.savefig(out_path_pdf + "Elements/alfa{}/".format(alfa) + filename + ".pdf")
                    plt.savefig(out_path_png + "Elements/alfa{}/".format(alfa) + filename + ".png", dpi = 600)


# Error plot
if conv_plot == True:
    for k in range(k_max):
        
        plt.figure(k +k_max*i_max*j_max +1)
        plt.xlabel("Element size")
        plt.ylabel(r"$||\epsilon_{\varepsilon}||_E$ relative error")
        # plt.ylabel("relative error")
        plt.grid(which = 'both')
        for i in range(i_max):
            plt.loglog(h_sizes, e[k,i,:], '-'+markers[i])
        plt.legend([r"$\alpha$ = " + "{}".format(alfa) for alfa in alfa_range])
        if mesh.elements[0,1] == 1:
            plt.title("Convergence analysis:\n{}-node isoparametric elements".format([2,3][k]))
            filename = "conv_{}-node_isoparemtric".format([2,3][k])
                    
        elif mesh.elements[0,1] == 2:
            plt.title("Convergence analysis:\n{}-node elements: closed formulation".format([2,3][k]))
            filename = "conv_{}-node_closedform".format([2,3][k])
            
        if save_conv == True:
            plt.savefig(out_path_eps + "Convergence/" + filename + ".eps")
            plt.savefig(out_path_pdf + "Convergence/" + filename + ".pdf")
            plt.savefig(out_path_png + "Convergence/" + filename + ".png")
            
        
plt.show()    
