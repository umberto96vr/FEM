import numpy as np

class Material():
    
    """ Contains information about the elements 
    material properties"""
    
    def __init__(self, k_spring = None):

        """k_spring is an array containing information about
        the elements stiffness"""
        self.k_spring = k_spring


