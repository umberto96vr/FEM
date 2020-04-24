# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 19:51:19 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys

import numpy as np
from numpy.linalg import inv

from   FEM.CustomErrors import FemError
import FEM.FEM_utilities as FEM_utilities
import FEM.material_data as material_data

# import time

# Additional utilities
integration = FEM_utilities.Integration()

#----------------------------------------------------------

def stiffness_matrix(e, mesh, MaterialSets, parameters):
    """
    Return the local stiffness matrix.
    
    Parameters
    ----------
    e:           element number
    mesh:        object containing mesh info
    material:    object containing material info
    utilities:   object containing integration quadrature scheme info
                 and shape functions derivatives
    
    Returns
    ---------
    k : local stiffness matrix
    
    """
    # Loading element type and material set number
    elType = mesh.elementType(e)
    elMatSet  = mesh.elements[e,1]
    
    # Loading stiffness matrix evaluation infos
    evaluation, domain, rule, points = material_data.stiffnes_matrix_eval_info(MaterialSets, elMatSet, elType)
    
    # Initializing variables
    strainSize = mesh.strain_sizes[mesh.d-1]
    D = np.zeros([strainSize, strainSize])
    
    
    if evaluation == 'closed form':  
                
        if elType == 'spring':
            
            k_s = material_data.elastic_properties(MaterialSets, elMatSet, elType)
    
            return k_s*np.array([(1, -1), (-1, 1)])
        
        elif elType == 'bar':
            
            # Loading cross section A and Young's modulus E
            A = material_data.geometric_properties(MaterialSets, elMatSet, elType)
            E = material_data.elastic_properties(MaterialSets, elMatSet, elType)
                        
            # Computing element length                    
            L = mesh.elem_coordinates(e)[1] - mesh.elem_coordinates(e)[0]
            
            if L < 0:
                
                print("Error: bar lenght is negative!")
                sys.exit()
                
            else:
                
                return ((E*A)/L)*np.array([(1, -1), (-1, 1)])                                      
            
        else:
            print("Error: element number {} not yet avaiable!".format(mesh.elements[e,0]))
            sys.exit()
            
    elif evaluation == 'numerical integration':
        
        if elType == 'bar':
            
            # Loading quadrature scheme weights and evaluation_points
            (w, int_p) = integration.quadrature_rule(rule, domain, points)
            
            # Loading elastic properties
            A = material_data.geometric_properties(MaterialSets, elMatSet, elType)
            E = material_data.elastic_properties(MaterialSets, elMatSet, elType)
            
            # Defining volume factor
            Vf = A
            
            # Defining the stiffness tensor D
            D[0,0] = E
            
            L = mesh.elem_coordinates(e)[1] - mesh.elem_coordinates(e)[0]
            
            # Initializing local stiffness matrix
            elementDofs = mesh.NodesElement * mesh.dofsNode
            k = np.zeros([elementDofs, elementDofs])
            
            if L < 0:
                
                print("Error: bar lenght is negative!")
                sys.exit()
            else:
                
                # Loop for numerical integration
                
                for i in range(len(int_p)):
                    
                    # Computing shape functions derivatives vector                     
                    (detJ, B) = FEM_utilities.shape_functions(mesh.elem_coordinates(e), int_p[i], 
                                                  mesh.NodesElement, mesh.d)

                    # i-step numerical contribution to stiffness matrix evaluation
                    # jacobian * dot(B_transposed, B) * Young modulus * Area * weight
                    
                    # dN has to be expressed in terms of the global coordinate system (x,y,z)
                    # The integrand function is F(x,y)*|J|*d_xi * d_eta that is
                    #   B(x,y)^T * D * B(x,y) * |J| * d_xi * d_eta 
                    # = B(xi,eta)^T * (J^-1)^T* D * B(xi,eta) * J^T * |J| * d_xi * d_eta 
                    
                    k += detJ * B.T @ D @ B * Vf * w[i]
                     
            return k  
            
        else:
            raise FemError("stiffness matrix evaluation key '{}' is not defined!".format(evaluation))


#----------------------------------------------------------

def DofMap(e, mesh):
    """
    Return a vector containing the global dof numbers associated\n
    with the local dofs for the given element\n
    
    Parameters
    ----------
    e : element number\n
    mesh: instance of Mesh() class containing mesh infos
    
    Returns
    ---------
    dof : dof map for the given element
    """   

    NodesInElement  = mesh.nodesInElement(e)
    NodesPerElement = mesh.NodesElement
       
    # dof vector generation via list comprehension
    dof = [NodesInElement[i] for i in range(NodesPerElement)]
    return dof


#----------------------------------------------------------

def assemble(K,k,dof):
    """
    Assemble in the global stiffness matrix "K" the local stiffness 
    matrix "k" associated with the given element associated to the dofs 
    map "dof"
    
    Parameters
    ----------
    K  : Global stiffness matrix \n
    k  : local stiffness matrix  \n
    dof: dof map for the given element  \n
        
    Returns
    -------
    K  : Assembled global stiffness matrix
    """
    # Computing the element number of degrees of freedom
    elementDofs = len(dof)
    
    if elementDofs <= 2:
        
        for i in range(elementDofs):
            ii = dof[i]
            for j in range(elementDofs):
                jj = dof[j]
                K[ii,jj] += k[i,j]
    else:
        
        K[np.ix_(dof,dof)] += k

    return K


