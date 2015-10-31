__author__ = 'Nasser'

import numpy as np
from scipy.sparse.linalg import spsolve
from scipy import sparse
import matplotlib.pyplot as plt
import gmsh as gmsh
import element2dof as element2dof
import assemble2dof as assemble2dof
import plotter as plotter
import boundaryconditions2dof as boundaryconditions2dof
import processing as processing


def solver(mesh_name, material, body_forces, traction_imposed,
           displacement_imposed,
           plot_undeformed, plot_stress, plot_deformed, print_report, **kwargs):

    mesh = gmsh.parse(mesh_name)

    ele = element2dof.Matrices(mesh)

    s = mesh.surfaces
    mat_dic = {s[i]: material[j] for i, j in enumerate(material)}

    ele.stiffness(mat_dic)

    ele.body_forces(body_forces)

    ele.mass(mat_dic)

    k = assemble2dof.globalMatrix(ele.K, mesh)

    m = assemble2dof.globalMatrix(ele.M, mesh)
    print(np.size(m), m)

    p0q = assemble2dof.globalVector(ele.P0q, mesh)

    p0t = boundaryconditions2dof.neumann(mesh, traction_imposed)

    p0 = p0q + p0t

    km, p0m = boundaryconditions2dof.dirichlet(k, p0, mesh,displacement_imposed)

    ks = sparse.csc_matrix(km)

    u = spsolve(ks, p0m)

    if print_report['U'] == True:
        processing.print_dof_values(u, mesh, 'Displacement')

    ele.nodal_forces(u)
    p_node = assemble2dof.globalVector(ele.pEle, mesh)

    sNode, sEle, eEle = processing.stress_recovery_simple2(mesh, u, mat_dic,
                                                           ele.e0)
    print(mesh.num_ele)
    if print_report['stress'] == True:
        processing.print_ele_values(sEle, mesh, 'Stress')

    if print_report['strain'] == True:
        processing.print_ele_values(eEle, mesh, 'Strain')

    principal_max = processing.principal_stress_max(sNode[0], sNode[1], sNode[2])
    principal_min= processing.principal_stress_min(sNode[0], sNode[1], sNode[2])

    dpi = 90
    magnification = plot_deformed['DeformationMagf']

    # PLOTTER CONTOUR MAPS
    if plot_stress['s11'] == True:
        plotter.tricontourf(sNode[0]/10**3, mesh,
                            'Stress 11 (kPa)', 'spring', dpi)

    if plot_stress['s22'] == True:
        plotter.tricontourf(sNode[1]/10**3, mesh,
                            'Stress 22 (kPa)', 'cool', dpi)

    if plot_stress['s12'] == True:
        plotter.tricontourf(sNode[2]/10**3, mesh,
                            'Stress 12 (kPa)', 'hsv', dpi)

    if plot_stress['sPmax'] == True:
        plotter.tricontourf(principal_max/10**3, mesh,
                            'Stress Principal Max (kPa)', 'autumn', dpi)

    if plot_stress['sPmin'] == True:
        plotter.tricontourf(principal_min/10**3, mesh,
                            'Stress Principal min (kPa)', 'winter', dpi)

    # PLOTTER DRAW UNDEFORMED SHAPE, ELEMENTS, LABELS, BC
    if plot_undeformed['Domain'] == True:
        plotter.draw_domain(mesh, 'Case Study', dpi, 'k')

    if plot_undeformed['Elements'] == True:
        plotter.draw_elements(mesh, 'Case Study', dpi, 'k')

    if plot_undeformed['ElementLabel'] == True:
        plotter.draw_elements_label(mesh, 'Case Study', dpi)

    if plot_undeformed['EdgesLabel'] == True:
        plotter.draw_edges_label(mesh, 'Case Study', dpi)

    if plot_undeformed['NodeLabel'] == True:
        plotter.draw_nodes_label(mesh, 'Case Study', dpi)

    if plot_undeformed['SurfaceLabel'] == True:
        plotter.draw_surface_label(mesh, 'Case Study', dpi)


    # PLOTTER DEFORMED SHAPE
    if plot_deformed['DomainUndeformed'] == True:
        plotter.draw_domain(mesh, 'Deformed Shape', dpi, 'SteelBlue')

    if plot_deformed['ElementsUndeformed'] == True:
        plotter.draw_elements(mesh, 'Deformed Shape', dpi, 'SteelBlue')

    if plot_deformed['DomainDeformed'] == True:
        plotter.draw_deformed_domain(mesh, u, 'Deformed Shape', dpi,
                                     magnification, 'Tomato')

    if plot_deformed['ElementsDeformed'] == True:
        plotter.draw_deformed_elements(mesh, u, 'Deformed Shape', dpi,
                                       magnification, 'Tomato', 1)

    plt.show()
