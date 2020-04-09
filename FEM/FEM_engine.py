import numpy as np
from numpy.linalg import inv
import time
#----------------------------------------------------------

def stiffness_matrix(k_s):
    """
    Return the local stiffness matrix for a truss
    element with the given spring stiffness
    
    Parameters
    ----------
    k_s: Spring stiffness
    
    Returns
    ---------
    K_spring : Spring local stiffness matrix
    """

    return k_s*np.array([(1, -1), (-1, 1)])

#----------------------------------------------------------

def DofMap(e, mesh):
    """
    Return a vector containing the global dof numbers associated
    with the local dofs for the given element
    
    Parameters
    ----------
    e : element number
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
    K  : Global stiffness matrix
    k  : local stiffness matrix
    dof: dof map for the given element
        
    Returns
    -------
    K  : Assembled global stiffness matrix

    """

    elementDofs = dof.size
    #print("dof vector dimension:", elementDofs)

    for i in range(elementDofs):
        ii = dof[i]
        for j in range(elementDofs):
            jj = dof[j]
            #print("local k[{},{}] goes into global K[{},{}]".format(i,j,ii,jj))
            #print("")
            K[ii,jj] += k[i,j]

    return K

#----------------------------------------------------------


def solve(K,F,ConstrainedDofs):
    """
    Solve the structural problem 
    F = K*U
    
    Parameters
    ----------
    K             : Global Stiffness
                    Matrix.
    F             : Load Vector.
    CostraineDofs : Vector containing
                    the constrained DoFs
                    numbers.

    Returns
    -------
    U             : Global DoFs vector
    R             : Global reactions vector

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