import numpy as np

#-----------------------------------------------------------------

class Mesh():
    """ Contains information about the discretization"""

    def __init__(self,ElementsNumber, Nodes, ...
                  NodesPerElement, dofsPerNode, ...
                  connect_table)
        
        self.Elements     = ElementsNumber
        self.Nodes        = Nodes
        self.NodesElement = NodesPerElement
        self.dofsNode     = dofsPerNode

        self.elements     = connect_table
