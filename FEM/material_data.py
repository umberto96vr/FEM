"""
Created on Fri Mar 20 14:39:22 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

#import numpy as np
import sys

class Material():
    
    """ 
    Contains information about the elements 
    material properties.
    
    materialSet is a nested dictionary with the following
    structure.
    {Material number: {Element type: ..., Material Properties: ... }
     
    Currently implemented material types are:
    -'1D Spring'--> Mat. Properties: 'Stiffness'
    -'1D Bar' ----> Mat. Properties: 'Area', 'Young Modulus'
      
    """
    def __init__(self, materialSet = None):
        
        self.materialSet = materialSet
        
    def parameters(self, elMat, elType):
        """       
        Parameters
        ----------
        elMat : Material number.
        elType : Element type.

        Returns
        -------
        Element parameters:
            1D Spring ---> spring stiffness
            1D Bar    ---> bar axial stiffness

        """
        
        if self.materialSet[str(elMat)]['Type'] == '1D Spring':
            
            return self.materialSet[str(elMat)]['Stiffness']
        
        elif self.materialSet[str(elMat)]['Type'] == '1D Bar':
            
            A = self.materialSet[str(elMat)]['Area']
            E = self.materialSet[str(elMat)]['Young Modulus']
            
            return E*A
        
        else:
            print('Error: wrong material type!')
            sys.exit()


