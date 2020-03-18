import numpy as np

#----------------------------------------------------------

def stiffness_matrix(springStiffness):
    """
    Return the local stiffness matrix for a truss
    element with the given spring stiffness
    """

    return springStiffness*np.array([(1, -1), (-1, 1)])

def DofMap(NodesInElement, nodesPerElement):
    """
    Return a vector containing the global dof numbers associated
    with the local dofs for the given element

    NodesInElement : nodes in the specified element
    NodesPerElement: number of nodes per element
    """   

    dof = np.zeros((nodesPerElement), dtype = np.int16)

    for i in range(nodesPerElement):
        dof[i] = NodesInElement[i]
    return dof

def assemble(K,k,dof):
    """
    Assemble in the global stiffness matrix "K" the local stiffness 
    matrix "k" associated with the given element associated to the dofs 
    map "dof"

    k:   local stiffness matrix
    dof: dof map for the given element
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