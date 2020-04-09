# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:53:10 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""


class FemError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class FemEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class MeshEngineError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)
        
class MaterialDataError(Exception):
    def __init__(self, text):
        self.text = text
        Exception.__init__(self, text)