# -*- coding: utf-8 -*-
"""
Created on Sat May 16 14:52:12 2020

@author: Umberto Bettinardi
         bettinardi96@gmail.com
"""
import sys
import time
sys.path.append("../../source")

import numpy as np
import meshio
from FEM import mesh_engine 
from FEM import BC_engine
import FEM.solvers.solver as solver
# from FEM.post.plotter_2d import plot_points
# import FEM.post.post2d as post2d
from assignment_support_functions import u_x, u_y, s_x, s_y, s_xy, change_mesh_size  
import matplotlib.pyplot as plt
import datetime

data_file = datetime.datetime.now().strftime("%Y-%m-%d_%H_%M")
figure_folder = "./figures/"
mesh_file = "kirsch_plate"

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#---------------------------------------------------------------

# Problem data
sigma_applied = 100
Young = 2.1*1e7
Poisson = 0.29
radius = 2.5
thikness = 1

#---------------------------------------------------------------
# MODEL PARAMETERS
#---------------------------------------------------------------

print("\n--> Pre-processing...\n")
print("Initializing model parameters...")

Procedures = {"solver": {"type": "linear"}, 
              
              "postprocessing":{"undeformed mesh": {"output": "screen"},
                                "deformed mesh": {"magnification factor": 50, "output": "screen"},
                                "axial stress": {"output": "screen"}
                                }
              }
#---------------------------------------------------------------

# Material info
print("Defining material sets...")

MaterialSets = {
    'piastra': {'element': 'triangle',
          'plane deformation'    : 'plane strain',
          'material behavior': 'isotropic linear elastic',
          'elastic properties'   : {"Young modulus" : Young, "poisson ratio" : Poisson},
          'geometric properties' : {'thikness': thikness},
          'stiffness matrix'     : {'evaluation': 'closed form'}
          
          }
    }
#---------------------------------------------------------------

H = [3.5, 1.7]#, 0.71, 0.35, 0.17, 0.14, 0.13]#, 0.12, 0.11, 0.1]
ERRORE  = []
ELEMEN  = []
DOFS    = []
REL_ERR = []

for h in H:

    print("\nMesh size: {}\n".format(h))
        
    # Mesh info
    print("Meshing...")
    
    # Upgrading mesh size in .geo file
    change_mesh_size(mesh_file, h)   # see assignment_support_function.change_mesh_size()
    mesh = mesh_engine.GMSH(mesh_file)  
    #--------------------------------------------------------------------------------------------------
    
    # BCs enforcement    
    print("Applying BCs...\n")
    BCs = BC_engine.BoundaryConditions()
    
    # Non-homogeneous Dirichlet BCs enforcement on the boundary
    name = "EnforcedDisplacement"
    tag = mesh.field_data[name][0]    
    on_boundary = np.nonzero(mesh.cell_data_dict["gmsh:physical"]["line"] == tag)[0]
    nodes_enf_disp = np.unique(mesh.cells_dict["line"][on_boundary])
    
    BCdata = []
 
    for i in nodes_enf_disp:
        
        x = mesh.points[i,0]
        y = mesh.points[i,1]
        
        ux = u_x(x, y, radius, Young, Poisson, sigma_applied)
        uy = u_y(x, y, radius, Young, Poisson, sigma_applied)        
                
        BCdata.append([int(i), "Dirichlet", 0, ux ])
        BCdata.append([int(i), "Dirichlet", 1, uy ])
        
    # In case symmetries are exploited corresponding BCs are enforced
    if mesh_file == 'kirsch_plate':
        BCdata += BCs.set_bc("Dirichlet", 'x-symmetry', mesh, [1], [0.])
        BCdata += BCs.set_bc("Dirichlet", 'y-symmetry', mesh, [0], [0.])
        
    BCs.data = BCdata
    
    #---------------------------------------------------------------
    # SOLVER
    
    print("--> Solving...\n")

    start = time.time()    
    U, R, K = solver.run(mesh, BCs, MaterialSets, Procedures)    
    end = time.time()    
    
    print("Total time: {}s".format(end-start))    
    print("\n-- Post-processing...\n")
        
    U_print = U.copy().reshape(len(mesh.points),mesh.d)
    print("   . solution U:\n{}\n".format(U_print))
    
    # -----------------------------------------------------------------------------
    # Stress computing 
    
    # Retrieving true displacements
    U_true = U_print.copy()
    for i in range(len(mesh.points)):
        # Points coordinates
        x_point = mesh.points[i,0]
        y_point = mesh.points[i,1]
        # Points displacements
        uu_x = u_x(x_point, y_point, radius, Young, Poisson, sigma_applied)
        uu_y = u_y(x_point, y_point, radius, Young, Poisson, sigma_applied)
        U_true[i,0] = uu_x
        U_true[i,1] = uu_y
    
    # Retrieving element number
    ElementNumber = len(mesh.elements) 
        
    # Retrieving strain size
    if mesh.d == 1:
        strainComponents = 1
    elif mesh.d == 2:
        strainComponents = 3
    elif mesh.d == 3:
        strainComponents = 6
        
    # Initializin results arrays
    stress      = np.zeros([ElementNumber, strainComponents])
    stress_true = np.zeros([ElementNumber, strainComponents])
    von_mises   = np.zeros(len(stress))
    elem_area   = np.zeros(len(stress))
    err_H0 = 0
    err_S  = 0
    
    for e in range(len(mesh.elements)):
        
        #dofs = int(mesh.dofsNode*len(mesh.elements[e]))
        
        # Loading element type and material set number
        elMatSet  = mesh.elementMaterialTag[e]
        key = MaterialSets[str(elMatSet)]
        
        # Loading thikness t, Young's modulus E and Poisson ration v.
        E = key['elastic properties'  ]['Young modulus']
        v = key['elastic properties'  ]['poisson ratio']
        t = key['geometric properties'  ]['thikness']
        
        
        n1 = mesh.elements[e][0]
        n2 = mesh.elements[e][1]
        n3 = mesh.elements[e][2]
        
        x1 = mesh.points[n1,0]
        y1 = mesh.points[n1,1]
        x2 = mesh.points[n2,0]
        y2 = mesh.points[n2,1]
        x3 = mesh.points[n3,0]
        y3 = mesh.points[n3,1]
        
        y23 = y2 - y3
        y31 = y3 - y1
        y12 = y1 - y2
        x32 = x3 - x2        
        x13 = x1 - x3
        x21 = x2 - x1
        
        A = 0.5*(x13*y23 - y31*x32)
        detJ = 2*A          
        elem_area[e] = A
        
        B = (1/detJ)*np.array([[y23,   0, y31,   0, y12,   0],
                               [  0, x32,   0, x13,   0, x21],
                               [x32, y23, x13, y31, x21, y12]])
        
        # Stiffness tensor evaluation
        if key['plane deformation'] == 'plane stress':
            
            D = (E/(1 - v**2))*np.array([[1, v,         0],
                                         [v, 1,         0],
                                         [0, 0, 0.5*(1-v)]])
        
        elif key['plane deformation'] == 'plane strain':
            
            D = (E/((1 - 2*v)*(1+v)))*np.array([[(1-v),     v,       0],
                                                [    v, (1-v),       0],
                                                [    0,     0, (1-2*v)/2.]])
            
        U_ = U_print[mesh.elements[e]].reshape([6,1])
        #U_ = U_true[mesh.elements[e]].reshape([6,1])
        
        sigma = np.linalg.multi_dot([D, B, U_])
        
        stress[e] = sigma.reshape(1,3)
        
        von_mises[e] = np.sqrt(sigma[0,0]**2 + sigma[1,0]**2 - sigma[0,0]*sigma[1,0] + 3*sigma[2,0]**2)
        
        # Computing element centroid
        x_c = (x1 + x2 + x3)/3
        y_c = (y1 + y2 + y3)/3
        
        # Computing "real" stresses
        stress_true[e,0] = s_x(x_c, y_c, radius, sigma_applied)
        stress_true[e,1] = s_y(x_c, y_c, radius, sigma_applied)
        stress_true[e,2] = s_xy(x_c, y_c, radius, sigma_applied)
        
        # Computin error
        n = 3
        def xx(xi,eta,X):
                
                return X[0] + (X[1]-X[0])*xi + (X[2]-X[0])*eta       
        
        def yy(xi,eta,Y):
            
            return Y[0] + (Y[1]-Y[0])*xi + (Y[2]-Y[0])*eta
        
        quad_scheme_GL = {
                            1: np.array([
                                        [1],
                                        [0.3333333333333333],
                                        [0.3333333333333333]]),
                            3: np.array([
                                         [0.3333333333333333,  
                                          0.3333333333333333,  
                                          0.3333333333333333],
                                         [0, 0.5, 0.5],
                                         [0.5, 0, 0.5]])}
        
        I = 0   # Initialize integral value
        Is = 0
        
        X = mesh.points[mesh.elements[e],0]
        Y = mesh.points[mesh.elements[e],1]
            
        # Computation loop
        
        for i in range(n):
            
            xi_i = quad_scheme_GL[n][1][i]
            eta_i = quad_scheme_GL[n][2][i]
            w_i = quad_scheme_GL[n][0][i]
            
            x_c = xx(xi_i, eta_i, X)
            y_c = yy(xi_i, eta_i, Y)
            
            stress_true[e,0] = s_x(x_c, y_c, radius, sigma_applied)
            stress_true[e,1] = s_y(x_c, y_c, radius, sigma_applied)
            stress_true[e,2] = s_xy(x_c, y_c, radius, sigma_applied)
            
            v = (stress_true[e,:] - stress[e,:]).reshape(3,1)
            vs = (stress_true[e,:]).reshape(3,1)
            
            I +=  0.5*w_i* float((v.T @ v))
            Is +=  0.5*w_i* float((vs.T @ vs))
        
        I *= detJ*t
        Is *= detJ*t
        err_H0 += I
        err_S  += Is
        
    err_H0 = np.sqrt(err_H0)
    err_S = np.sqrt(err_S)
    REL_ERR.append(err_H0/err_S)
    ERRORE.append(err_H0)
    ELEMEN.append(len(mesh.elements))
    DOFS.append(len(K))
          
    F = np.matmul(K,U).reshape(len(mesh.points),mesh.d)
    #--------------------------------------------------------------------------------------------------
    
    # PARAVIEW .vtk file generation
    
    # Generate additional data for meshio for .vtk file generation
    mesh.point_data = {'Displacements_FEM': U.reshape(len(mesh.points),mesh.d), 
                       'Displacements_true': U_true,
                       "Forces" : F}   
    mesh.cell_data = {"stress_true": stress_true,"stress_fem": stress, "von_mises": von_mises}      
    cells = {'triangle': mesh.elements}
    
    # Write the vtk file
    meshio.write_points_cells("./fileVtk/{}/".format(
        MaterialSets['piastra']['plane deformation'].replace(" ","_"))+ mesh_file +"_e{}.vtk".format(
            len(mesh.elements)), mesh.points, cells, mesh.point_data, mesh.cell_data, binary=False)
    
    #-------------------------------------------------------------------------------------------------- 
    
    #Error norm computing
    E = 0
    for e in range(len(mesh.elements)):
    
        v = (stress_true[e,:] - stress[e,:]).reshape(3,1)
        
        # Integral - 2nd order Gaussian quadrature
    
        E += elem_area[e] * float((v.T @ v))
    
    err_H0_2 = np.sqrt(E)
    
    E = 0
    for i in range(len(mesh.points)):
        
        v = (U_print[i,:] - U_true[i,:]).reshape(2,1)
        E += float((v.T @ v))
    err_U = np.sqrt(E)
    
    print("\nerr_H0 : {}".format(err_H0))
    
conv_rate = -np.log(ERRORE[-1]/ERRORE[-2])/np.log(DOFS[-1]/DOFS[-2])
#--------------------------------------------------------------------------------------------------

if len(ELEMEN)>=2:
    #--------------------------------------------------------------------------------------------------        
    conv_rate = -np.log(ERRORE[-1]/ERRORE[-2])/np.log(DOFS[-1]/DOFS[-2])
    #--------------------------------------------------------------------------------------------------    
    
    # Writing results to file
    result_file = open('./results/'+mesh_file+'_results_'+data_file+'.txt', 'w')
    
    result_file.write("{}.geo - results\n\n".format(mesh_file))
    result_file.write("@author: Umberto Bettinardi\n"+
                      "         bettinardi96@gmail.com\n\n")
    result_file.write("Analysis date: {}\n\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M")))
    result_file.write("\t{:.8s}\t  {:.8s}\t    {:.8s}\t\t{:.8s}\t {:.9s}\n".format("h-size","Elements", "Dofs", "err_H0", "rel_err"))
    for i in range(len(ELEMEN)):
        result_file.write("{:8.3f}\t{:8d}\t{:8d}\t{:8.6f}\t{:9.6f}\n".format(H[i],ELEMEN[i],DOFS[i],ERRORE[i], REL_ERR[i]))
    result_file.write("\nConvergence rate:\t{}".format(conv_rate))    
    result_file.close()    
    #-------------------------------------------------------------------------------------------------- 
    
    # Plotting convergence diagramm
    fig, ax = plt.subplots()
    ax.grid(which = 'both')
    ax.set_title("Convergence analysis: T3 elements.")
    ax.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
    ax.set_xlabel("Number of dofs")
    ax.loglog(ELEMEN,ERRORE,'-*b')   
    textstr = "Convergence rate: {:4.3f}".format(conv_rate)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.5, 0.9, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    fig.savefig(figure_folder+"convergence_{}.eps".format(mesh_file), dpi = 600)
    fig.savefig(figure_folder+"convergence_{}.png".format(mesh_file), dpi = 600)
    fig.savefig(figure_folder+"convergence_{}.pdf".format(mesh_file))
    
    #--------------------------------------------------------------------------------------------------        
    conv_rate = np.log(ERRORE[-1]/ERRORE[-2])/np.log(H[-1]/H[-2])
    #--------------------------------------------------------------------------------------------------  
    fig2, ax2 = plt.subplots()
    ax2.grid(which = 'both')
    ax2.set_title("Convergence analysis: T3 elements.")
    ax2.set_ylabel(r"$ ||\sigma_{e} - \sigma_{h}||_{H^0}$")
    ax2.set_xlabel("Element size")
    ax2.loglog(ELEMEN,ERRORE,'-*b')
    textstr = "Convergence rate: {:4.3f}".format(conv_rate)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax2.text(0.5, 0.9, textstr, transform=ax2.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
    fig2.savefig(figure_folder+"convergence_{}_H.eps".format(mesh_file), dpi = 600)
    fig2.savefig(figure_folder+"convergence_{}_H.png".format(mesh_file), dpi = 600)
    fig2.savefig(figure_folder+"convergence_{}_H.pdf".format(mesh_file))
    
    
    plt.show()

#--------------------------------------------------------------------------------------------------