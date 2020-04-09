"""
Created on Fri Mar 20 14:45:31 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
#import numpy as np

#-----------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization"""

    def __init__(self,ElementsNumber = None, Nodes = None, NodesPerElement = None,
                 dofsPerNode = None, connect_table = None):

        """
        Return a vector containing the global dof numbers associated
        with the local dofs for the given element
        
        Attributes
        ----------
        ElementsNumber : number of elements
        Nodes: number of nodes per element
        NodesPerElement
        dofsPerNode
        conncet_table: requires an array of dimension n x 4 where n 
                       is the number of elements. 
                       Each row has the following structure: 
                       element type | material | 1st node in global ref | 2nd node in global ref
        """
        
        self.Elements     = ElementsNumber
        self.Nodes        = Nodes
        self.NodesElement = NodesPerElement
        self.dofsNode     = dofsPerNode
        self.elements     = connect_table 
      
    
    def dofs(self):
        """
        Return the total number of degrees of freedom

        """
        
        return self.Nodes*self.dofsNode
    
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