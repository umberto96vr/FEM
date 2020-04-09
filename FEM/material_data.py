"""
Created on Fri Mar 20 14:39:22 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

#import numpy as np
# import sys
from FEM.CustomErrors import FemError, MaterialDataError

class Material():
    
    """ 
    Contains information about the elements 
    material properties.
    
    materialSet is a nested dictionary with the following
    structure:\n
    {Material number: {Element type: ..., Elastic Properties: {...}, Geometric Properties: {...}}
     
    Currently implemented material types are:\n
    -'spring'--> Elast. Properties: 'stiffness'\n
    -'bar' ----> Elast. Properties: 'Young Modulus'\n
                 Geome. Properties: 'area'
      
    """
    
    #------------------------------------------------------------------------------------------
    
    def __init__(self, materialSet = None):
        
        self.materialSet = materialSet
     
    #------------------------------------------------------------------------------------------
       
    def elastic_properties(self, elMatSet, elType):
        """       
        Return the elastic properties\n
        Parameters
        ----------
        elMatSet : Material set number.\n
        elType   : Element type.\n

        Returns:
        -------
        Element parameters:\n\t
            spring ---> spring stiffness\n\t
            bar    ---> Young modulus\n

        """

        key = self.materialSet[str(elMatSet)]

        
        if key['element'] == 'spring' and elType == 'spring':
            
            return key['elastic properties']['stiffness']
        
        elif key['element'] == 'bar' and elType == 'bar':
            
            return key['elastic properties']["Young modulus"]
        
        else:
            
            error_msg = 'Error: either {} is not a valid element type or\n)'.format(elType)
            error_msg += '       {} is not a valid material number!\n'.format(elMatSet)
            
            raise MaterialDataError(error_msg)
            
            # print('Error: either {} is not a valid element type or\n'.format(elType))
            # print('       {} is not a valid material number!\n'.format(elMatSet))
            # sys.exit()
    
    #------------------------------------------------------------------------------------------
            
    def geometric_properties(self, elMatSet, elType):
        """       
        Return the geometric properties\n
        Parameters
        ----------
        elMatSet : Material set number.\n
        elType   : Element type.\n

        Returns
        -------
        Element parameters:\n\t
            bar  ---> bar cross section

        """
        
        key = self.materialSet[str(elMatSet)]
        
        
        if key['element'] == 'bar' and elType == 'bar':
            
            return key['geometric properties']['area']
 
        else:
            
            error_msg = 'Error: either {} is not a valid element type or\n)'.format(elType)
            error_msg += '       {} is not a valid material number!\n'.format(elMatSet)
            
            raise MaterialDataError(error_msg)
            
            # print('Error: either {} is not a valid element type or\n'.format(elType))
            # print('       {} is not a valid material number!\n'.format(elMatSet))
            # sys.exit()
    
    #------------------------------------------------------------------------------------------
    
    def stiffnes_matrix_eval_info(self, elMatSet, elType):
        """
        Returns info concerning the evaluation of the
        stiffness matrix\n
        Parameters
        ----------
        elMatSet : Material set number.\n
        elType   : Element type. \n

        Returns
        -------
        evaluation : evaluation mode of the stiffness matrix. \n
        domain     : domain of integration. \n
        rule       : quadrature rule of integration. \n
        points     : number of integration points. \n

        """
        key = self.materialSet[str(elMatSet)]['stiffness matrix']
        
        if key['evaluation'] == 'closed form':
            
            evaluation = 'closed form'
            domain     = None
            rule       = None
            points     = None
 
        elif key['evaluation'] == 'numerical integration':
            
            evaluation = key['evaluation']
            domain     = key['domain']
            rule       = key['rule']
            points     = key['points']
        else:
            raise MaterialDataError('Error: {} is not a valid evaluation option!'.format(key['evaluation']))
            
            # print('Error: {} is not a valid evaluation option!'.format(key['evaluation']))
            # sys.exit()
        
        return evaluation, domain, rule, points