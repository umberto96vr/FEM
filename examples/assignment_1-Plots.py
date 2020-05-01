# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 08:54:21 2020

@author: Bettinardi Umberto
@mail:   umberto.bettinardi@studenti.unipd.it
         bettinardi96@gmail.com
"""
import numpy as np
import matplotlib.pyplot as plt

# System parameters
L = 10
P = 1
E = 1
A = 1
alfa = [0.5, 1, 2, 4]
k = [E*A*i**2 for i in alfa]

# Analysis domain
x = np.linspace(0,L,100)
y = np.zeros(shape = (len(x),len(alfa)))

for i in range(len(alfa)):
    y[:,i] = -P/(E*A*alfa[i]) * np.exp(-alfa[i]*x)
    
# Plotting

colors = ['r','b','g','m']
markers = ['*', 'o', 'x', '^']

plt.figure(1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif') 
for i in range(len(alfa)):
    style = '-'+colors[i]#+markers[i]
    
    plt.plot(x, y[:,i], style)
    
   
plt.xlabel(r'$x$', fontsize=12)
plt.ylabel(r'$u (x)$', fontsize=12)
plt.legend([r"$\alpha$ = {}".format(str(i)) for i in alfa])
plt.title(r"$Analytical\; solution$", fontsize=15)
plt.grid()
plt.savefig('.\\Plots\\analytical_sol.eps', dpi = 600)