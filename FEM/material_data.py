"""
Created on Fri Mar 20 14:39:22 2020
Modified on Fri Apr 20 09:47:51 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""

from FEM.CustomErrors import FemError, MaterialDataError

def elastic_properties(MaterialSet, elMatSet, elType):
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

    key = MaterialSet[str(elMatSet)]

    
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
        
def geometric_properties(MaterialSet, elMatSet, elType):
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
    
    key = MaterialSet[str(elMatSet)]
    
    
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

def stiffnes_matrix_eval_info(MaterialSet, elMatSet, elType):
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
    key = MaterialSet[str(elMatSet)]['stiffness matrix']
    
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