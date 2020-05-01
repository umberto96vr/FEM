# -*- coding: utf-8 -*-
"""
Created on Fri May  1 10:26:12 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np

v = np.random.randint(-5,5,50).reshape(-1,1)
A = v @ v.T

K = 0.5*(A + A.T)

K[:,[3, 5, 8]] = 0
K[[3, 5, 8], :] = 0
K[[3, 5, 8],[3, 5, 8]] = 1

print(np.linalg.det(K))
