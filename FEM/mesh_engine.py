"""
Created on Fri Mar 20 14:45:31 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
#import numpy as np
# import sys
from FEM.CustomErrors import FemError, MeshEngineError

#---------------------------------------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization\n
        Attributes: \n
        ----------
        ElementsNumber : number of elements. \n
        Nodes: number of nodes per element. \n
        NodesPerElement. \n
        dofsPerNode. \n
        connect_table: \n
        \t\t\t         requires an array of dimension n x 4 where n \n 
        \t\t\t         is the number of elements. \n
        \t\t\t         Each row has the following structure: \n
        \t\t\t         element type | material | 1st node in global ref | 2nd node in global ref
        """
        
    #------------------------------------------------------------------------------------------

    def __init__(self,ElementsNumber = None, Nodes = None, NodesPerElement = None,
                 dofsPerNode = None, connect_table = None, coordinates = None):
        
        self.Elements     = ElementsNumber
        self.Nodes        = Nodes
        self.NodesElement = NodesPerElement
        self.dofsNode     = dofsPerNode
        self.elements     = connect_table 
        self.coordinates  = coordinates
      
        
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
        for the specified element.

        """
        
        return self.elements[e][2:]
        
    #------------------------------------------------------------------------------------------
    
    def elementType(self,e):
        """
        Returns the element type for element e

        Parameters
        ----------
        e : element number\n
        
        Implemented elements:\n\t
            
            1   spring element\n\t
            2   bar element\n

        Returns
        -------
        element type

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