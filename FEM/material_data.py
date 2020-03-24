"""
Created on Fri Mar 20 14:39:22 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

#import numpy as np

class Material():
    
    """ Contains information about the elements 
    material properties"""
    
    def __init__(self, k_spring = None):

        """k_spring is an array containing information about
        the elements stiffness"""
        self.k_spring = k_spring


