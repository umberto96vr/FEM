import numpy as np

#-----------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization"""

    def __init__(self,ElementsNumber = None, Nodes = None, NodesPerElement = None,
                 dofsPerNode = None, connect_table = None):
        
        self.Elements     = ElementsNumber
        self.Nodes        = Nodes
        self.NodesElement = NodesPerElement
        self.dofsNode     = dofsPerNode

        self.elements     = connect_table
