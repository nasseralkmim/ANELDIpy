import eldypy.plotter as plotter
import matplotlib.pyplot as plt


def data(plot_stress, plot_deformed, plot_undeformed, mesh, s_ele, e_ele,
         s_node, principal_max, principal_min, u):

    dpi = 90
    magnification = plot_deformed['DeformationMagf']

    # PLOTTER CONTOUR MAPS
    if plot_stress['s11'] is True:
        plotter.tricontourf(s_node[0]/10**3, mesh,
                            'Stress 11 (kPa)', 'spring', dpi)

    if plot_stress['s22'] is True:
        plotter.tricontourf(s_node[1]/10**3, mesh,
                            'Stress 22 (kPa)', 'cool', dpi)

    if plot_stress['s12'] is True:
        plotter.tricontourf(s_node[2]/10**3, mesh,
                            'Stress 12 (kPa)', 'hsv', dpi)

    if plot_stress['sPmax'] is True:
        plotter.tricontourf(principal_max/10**3, mesh,
                            'Stress Principal Max (kPa)', 'autumn', dpi)

    if plot_stress['sPmin'] is True:
        plotter.tricontourf(principal_min/10**3, mesh,
                            'Stress Principal min (kPa)', 'winter', dpi)

    # PLOTTER DRAW UNDEFORMED SHAPE, ELEMENTS, LABELS, BC
    if plot_undeformed['Domain'] is True:
        plotter.draw_domain(mesh, 'Case Study', dpi, 'k')

    if plot_undeformed['Elements'] is True:
        plotter.draw_elements(mesh, 'Case Study', dpi, 'k')

    if plot_undeformed['ElementLabel'] is True:
        plotter.draw_elements_label(mesh, 'Case Study', dpi)

    if plot_undeformed['EdgesLabel'] is True:
        plotter.draw_edges_label(mesh, 'Case Study', dpi)

    if plot_undeformed['NodeLabel'] is True:
        plotter.draw_nodes_label(mesh, 'Case Study', dpi)

    if plot_undeformed['SurfaceLabel'] is True:
        plotter.draw_surface_label(mesh, 'Case Study', dpi)

    # PLOTTER DEFORMED SHAPE
    if plot_deformed['DomainUndeformed'] is True:
        plotter.draw_domain(mesh, 'Deformed Shape', dpi, 'SteelBlue')

    if plot_deformed['ElementsUndeformed'] is True:
        plotter.draw_elements(mesh, 'Deformed Shape', dpi, 'SteelBlue')

    if plot_deformed['DomainDeformed'] is True:
        plotter.draw_deformed_domain(mesh, u, 'Deformed Shape', dpi,
                                     magnification, 'Tomato')

    if plot_deformed['ElementsDeformed'] is True:
        plotter.draw_deformed_elements(mesh, u, 'Deformed Shape', dpi,
                                       magnification, 'Tomato', 1)

    plt.show()
