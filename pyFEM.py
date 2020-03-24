# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 17:38:10 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import numpy as np

from FEM import FEM_engine as engine
from FEM import mesh_engine as mesher
from FEM import material_data as material
from FEM import boundaries 

#------------------------------------------------------

mesh_info = mesher.Mesh()
material_info = material.Material()
load_info = boundaries.BoundaryConditions()

#------------------------------------------------------

# Mesh info
mesh_info.Elements = 4
mesh_info.Nodes = 4
mesh_info.NodesElement = 2
mesh_info.dofsNode = 1
mesh_info.elements = np.array([[0,1],[1,2],[1,2],[2,3]])

#------------------------------------------------------
# Material info
material_info.k_spring = np.array([20, 30, 2, 10])

#------------------------------------------------------
load_info.loads = np.array([0, 1, 0, -3])
load_info.constraints = np.array([0])

#------------------------------------------------------
#                   ASSEMBLY
#------------------------------------------------------

print("-------------------------------------------------------")
print("                     ASSEMBLY")
print("-------------------------------------------------------")

dofs = mesh_info.dofsNode*mesh_info.Nodes
K = np.zeros(shape = (dofs, dofs))

for e in range(mesh_info.Elements):
    
    k = engine.stiffness_matrix(material_info.k_spring[e])
    dof = engine.DofMap(mesh_info.elements[e], mesh_info.NodesElement)
    
    engine.assemble(K, k, dof)
print(K)
    
#------------------------------------------------------
#                   SOLVING
#------------------------------------------------------

print("-------------------------------------------------------")
print("                     SOLVING")
print("-------------------------------------------------------")

U, R = engine.solve(K,load_info.loads, load_info.constraints)

print(U)
print(R)
