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
        ElementsNumbert : nodes in the specified element
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