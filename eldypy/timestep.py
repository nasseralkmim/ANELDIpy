from eldypy import explicitfd
from eldypy import plotter
from eldypy import processing
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np


def iterations(period, steps, u0, v0, mesh, mat_dic, km, mm, p0m, interval):
    """ Perform all the steps defined for the period

    """
    dt = period/steps
    u1 = u0 + dt*v0
    u0_v = np.zeros(len(p0m))
    u1_v = np.zeros(len(p0m))
    u0_v[:] = u0
    u1_v[:] = u1

    G = nx.Graph()
    u_fig = plt.figure('Deformed Shape')
    #smax_fig = plt.figure('Maximum Principal Stress')
    #smin_fig = plt.figure('Minimum Principal Stress')

    u_ims = []
    smax_ims = []
    smin_ims = []

    u_initial_frame = plotter.draw_deformed_elements2(mesh, u0_v, 't', 200, 50, 'Tomato', G)
    u_ims.append((u_initial_frame,))


    s_node, s_ele, e_ele = processing.stress_recovery(mesh, u0_v, mat_dic)

    smax = processing.principal_stress_max(s_node[0], s_node[1],
                                                    s_node[2])
    smin = processing.principal_stress_min(s_node[0], s_node[1],
                                                    s_node[2])

    #smax_initial_frame = plotter.tricontourf2(smax/10**3, mesh, 'Stress Principal Max (kPa)', 'autumn', 50)

    #smin_initial_frame = plotter.tricontourf2(smin/10**3, mesh, 'Stress Principal min (kPa)', 'winter', 50)

    #smax_ims.append((smax_initial_frame,))
    #smin_ims.append((smin_initial_frame,))

    for t in range(steps):
        u_updt = explicitfd.scheme(u0_v, u1_v, mm, p0m, km, dt)
        u0_v = u1_v
        u1_v = u_updt

        s_node, s_ele, e_ele = processing.stress_recovery(mesh, u_updt, mat_dic)

        smax = processing.principal_stress_max(s_node[0], s_node[1],
                                                    s_node[2])
        smin = processing.principal_stress_min(s_node[0], s_node[1],
                                                    s_node[2])
        #smax_frame = plotter.tricontourf2(smax/10**3, mesh, 'Stress Principal Max (kPa)', 'autumn', 50)

        #smin_frame = plotter.tricontourf2(smin/10**3, mesh, 'Stress Principal min (kPa)', 'winter', 50)

        #smax_ims.append((smax_frame,))
        #smin_ims.append((smin_frame,))

        u_frame = plotter.draw_deformed_elements2(mesh, u_updt, 'Deformed Shape', 200, 50, 'Tomato', G)
        u_ims.append((u_frame,))

    u_ani = animation.ArtistAnimation(u_fig, u_ims, blit=True, interval=interval)

    #smax_ani = animation.ArtistAnimation(smax_fig, smax_ims, blit=True, interval=interval)
    #smin_ani = animation.ArtistAnimation(smin_fig, smin_ims, blit=True, interval=interval)
    #ani.save('smax_ani.mp4', dpi=300)
    plt.show()
