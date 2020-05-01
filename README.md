# SpritzFEM
## FEM program for linear-static analysis

This project is part of the course of COMPUTATIONAL MECHANICS OF MATERIAL for the MSc in Mechanical Engineering.

Goal is to build an all-around Fem code dealing with linear-elasticity (and eventually non-linear analysis and damage simulation).

The main engine is built around isoparametric elements that are integrated via Gauss-Legendre quadrature, 1D and 2D solver are under development, 3D will be eventually added in the future.

The code of the main engine is built on Python and resolve to other commercial open-source software for the meshing part and post-processing visualisation of results suchs as Paraview and gmsh
