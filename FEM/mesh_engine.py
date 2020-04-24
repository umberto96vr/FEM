"""
Created on Fri Mar 20 14:45:31 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
# import sys

import numpy as np

from FEM.CustomErrors import FemError, MeshEngineError

#---------------------------------------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization\n
        Parameters: 
        ----------
        ElementsNumber : number of elements. \n
        Nodes : number of nodes. \n
        NodesPerElement : number of nodes per element. \n
        dofsPerNode : number of degrees of freedom per node \n
        connect_table: \n
        \t\t\t         requires an array of dimension n x 4 where\n   
        \t\t\t         n is the number of elements. \n
        \t\t\t         Each row has the following structure: \n
        \t\t\t         element type | material | 1st node in global ref | 2nd node in global ref\n
        coordinates : element nodes coordinates.\n
        d : dimension of the problem.
        """
        
    #------------------------------------------------------------------------------------------

    def __init__(self):
        
        self.Elements     = None 
        self.Nodes        = None
        self.NodesElement = None
        self.dofsNode     = None
        self.elements     = None 
        self.coordinates  = None
        self.d            = None
        self.strain_sizes      = [1, 3, 6]
      
        
    #------------------------------------------------------------------------------------------
    
    def dofs(self):
        """
        Return the total number of degrees of freedom

        """
        
        return self.Nodes*self.dofsNode
        
    #------------------------------------------------------------------------------------------
    
    def nodesInElement(self,e):
        """
        Parameters
        ----------
        e : Element number.

        Returns
        -------
        Returns the global dofs mapping
        for the specified element e.

        """
        
        return self.elements[e,2:]
        
    #------------------------------------------------------------------------------------------
    
    def elementType(self,e):
        """
        Returns the element type for element e

        Parameters
        ----------
        e : element number\n
        
        Implemented elements:\n\t
            
            1   spring element \n\t
            2   bar element 

        Returns
        -------
        element type:\n\t
            
            1 ----> 'spring'\n
            2 ----> 'bar'
        

        """
        if self.elements[e,0] == 1:
            return "spring"
        
        elif self.elements[e,0] == 2:
            return "bar"
        else:
            raise MeshEngineError("Error: element number {} not yet avaiable!".format(self.elements[e,0]))
            # print("Error: element number {} not yet avaiable!".format(self.elements[e,0]))
            # sys.exit()
        
    #------------------------------------------------------------------------------------------
    
    def elem_coordinates(self, e):
        """
        Returns the nodes element coordinates xi, xj, for a bar element

        Parameters
        ----------
        e : element number.

        Returns
        -------
        X : coordinates vector.

        """
        
        return self.coordinates[self.nodesInElement(e)]
    
    #------------------------------------------------------------------------------------------
       
    def __str__(self):
        
        str_info = "------------------------------------------\n"
        str_info += "             MESH INFO\n"
        str_info += "------------------------------------------\n"
        str_info += "Number of elements:          {}\n".format(self.Elements)
        str_info += "Number of nodes:             {}\n".format(self.Nodes)
        str_info += "Number of nodes per elment:  {}\n".format(self.NodesElement)
        str_info += "Number of dofs per node:     {}\n\n".format(self.dofsNode)
        str_info += "Connectivity table:\n\n"
        str_info += "| Material |Element type  |   Node i    |Node j\n"
        for i in range(self.Elements):
            
            str_info += "|   {}      |        {}     |     {}       |     {}\n".format(
                self.elements[i,0],self.elements[i,1],self.elements[i,2],self.elements[i,3])
        return str_info