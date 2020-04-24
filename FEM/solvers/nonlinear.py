#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:12:21 2020

@author: angelo.simone@unipd.it
"""

class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

# -----------------------------------------------------------------------------

def nonlinear(mesh, BCs, MaterialSets, parameters):
    
    raise SolverError("Nonlinear solver not yet implemented!")