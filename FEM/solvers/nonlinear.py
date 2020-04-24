#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 11:09:56 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""

class SolverError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)

# -----------------------------------------------------------------------------

def nonlinear(mesh, BCs, MaterialSets, parameters):
    
    raise SolverError("Nonlinear solver not yet implemented!")