# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 18:30:32 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

# This reimplements gmsh/demos/boolean/boolean.geo in Python.
import meshio
import os
import numpy as np
import matplotlib.pyplot as plt

mesh_file = "bimaterial2"

os.system("gmsh -2 " + mesh_file + ".geo" + " -o " + mesh_file + ".msh")

mesh = meshio.read(mesh_file + ".msh")

triangle_cells = mesh.cells_dict["triangle"]
points = mesh.points
triangle_data = mesh.cell_data_dict["gmsh:physical"]["triangle"]

first_kind = len(triangle_data[triangle_data == 112])

plt.figure(1)
for i in triangle_cells[:first_kind+1]:
    X = points[[i],0][0]
    Y = points[[i],1][0]
    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    plt.plot(X,Y, '-r', linewidth = 0.5)
for i in triangle_cells[first_kind:]:
    X = points[[i],0][0]
    Y = points[[i],1][0]
    X = np.append(X, X[0])
    Y = np.append(Y, Y[0])
    plt.plot(X,Y, '-b', linewidth = 0.5)
    
size = triangle_cells.shape[0]

dirichlet = mesh.cell_sets_dict["Dirichlet"]["line"]

elem = np.zeros((size, 5), dtype = int)

for i in range(size):
    elem[i,:2] = np.array([triangle_data[i],6])
    elem[i,2:] = triangle_cells[i,:]

