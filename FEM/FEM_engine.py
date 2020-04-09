import numpy as np
from numpy.linalg import inv
from FEM.CustomErrors import FemError
import sys
# import time


class HostNotFound(Exception):
    def __init__( self, host ):
        self.host = host
        Exception.__init__(self, 'Host Not Found exception: missing %s' % host)


#----------------------------------------------------------

def stiffness_matrix(e, mesh, material):
    """
    Return the local stiffness matrix for a truss
    element with the given spring stiffness
    
    Parameters
    ----------
    e:          element number
    mesh:       object containing mesh info
    material:   object containing material info
    
    Returns
    ---------
    K_spring : Spring local stiffness matrix
    
    """
    
    e_type = mesh.elementType(e)
    elMatSet  = mesh.elements[e,1]
    evaluation, domain, rule, points = material.stiffnes_matrix_eval_info(elMatSet, e_type)
    
    if evaluation == 'closed form':    
        if e_type == 'spring':
            
            k_s = material.elastic_properties(elMatSet, e_type)
    
            return k_s*np.array([(1, -1), (-1, 1)])
        
        elif e_type == 'bar':
            
            A = material.geometric_properties(elMatSet, e_type)
            E = material.elastic_properties(elMatSet, e_type)
            
            n_i = mesh.nodesInElement(e)[0]
            n_j = mesh.nodesInElement(e)[1]
            
            L = mesh.coordinates[n_j] - mesh.coordinates[n_i]
            
            if L < 0:
                
                print("Error: bar lenght is negative!")
                sys.exit()
            else:
                return ((E*A)/L)*np.array([(1, -1), (-1, 1)])
                
            
            
            
        else:
            print("Error: element number {} not yet avaiable!".format(mesh.elements[e,0]))
            sys.exit()
            
    elif evaluation == 'numerical integration':
        raise FemError("stiffness matrix evaluation key '{}' hasn't been yet implemented!".format(evaluation))
        # print("Error: stiffness matrix evaluation key '{}' hasn't been yet implemented!".format(evaluation))
        sys.exit()

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
    
    

    dof = np.zeros((NodesPerElement), dtype = np.uint8)

    for i in range(NodesPerElement):
        
        if NodesInElement[i] > 255:
            
            print('Overflow error!')
        
        dof[i] = NodesInElement[i]
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

    elementDofs = dof.size

    for i in range(elementDofs):
        ii = dof[i]
        for j in range(elementDofs):
            jj = dof[j]
            K[ii,jj] += k[i,j]

    return K


#------------------------------------------------------------------------------------------
    
def solve(K,F,ConstrainedDofs):
    """
    Solve the structural problem  \n
    F = K*U
    
    Parameters
    ----------
    K             : Global Stiffness Matrix. \n 
    F             : Load Vector. \n
    CostraineDofs : Vector containing the constrained DoFs numbers. \n

    Returns
    -------
    U             : Global DoFs vector \n
    R             : Global reactions vector \n

    """
    
    
    if len(F) == len(K[:,0]):
        K_r = K.copy()
        
        for e in ConstrainedDofs:
            
            K_r[:,e] = 0
            K_r[e,:] = 0
            K_r[e,e] = 1
        U = np.matmul(inv(K_r),F)
        R = np.matmul(K,U)
        
        return U, R
    else:
        print('Error: K is a {}x{} matrix while \nF is {}x1! Shape mismatch!'.format(len(K[:,0]),len(K[:,0]),len(F)))
        return np.zeros((len(K[:,0]),1)), np.zeros((len(K[:,0]),1))