import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation


def draw_elements_label(mesh, name, dpi):
    plt.figure(name, dpi=dpi)

    for e in range(len(mesh.ele_conn)):
        x_element = (mesh.nodes_coord[mesh.ele_conn[e, 0], 0] +
                  mesh.nodes_coord[mesh.ele_conn[e, 1], 0] +
                          mesh.nodes_coord[mesh.ele_conn[e, 2], 0] +
                          mesh.nodes_coord[mesh.ele_conn[e, 3], 0])/4.

        y_element = (mesh.nodes_coord[mesh.ele_conn[e, 0], 1] +
                          mesh.nodes_coord[mesh.ele_conn[e, 1], 1] +
                          mesh.nodes_coord[mesh.ele_conn[e, 2], 1] +
                          mesh.nodes_coord[mesh.ele_conn[e, 3], 1])/4.


        plt.annotate(str(e), (x_element, y_element), size=9,
                     color='r')

        plt.axes().set_aspect('equal')
        plt.axes().autoscale_view(True, True, True)
        plt.margins(y=0.1, x=0.1, tight=False)
        plt.draw()


def draw_surface_label(mesh, name, dpi):

    c = mesh.nodes_coord
    plt.figure(name, dpi=dpi, frameon = False)

    X, Y = c[:, 0], c[:, 1]

    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    # adjust the node numbering order of an element
    if mesh.gmsh == 1.0:
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1

    i = 0
    for surface, lpTag in mesh.physicalSurface.items():
        xm = 0.0
        ym = 0.0
        for node in mesh.lineLoop[lpTag]:
            G2.add_edge(mesh.line[node][0],
                        mesh.line[node][1])
            xm += (mesh.nodes_coord[mesh.line[node][0], 0] +
                   mesh.nodes_coord[mesh.line[node][1], 0])
            ym += (mesh.nodes_coord[mesh.line[node][0], 1] +
                   mesh.nodes_coord[mesh.line[node][1], 1])

        xs, ys = xm/(2*len(mesh.lineLoop[lpTag])), ym/(2*len(mesh.lineLoop[
                                                                 lpTag]))
        plt.annotate(str(i), (xs, ys), size=9, color='g')
        i += 1

    positions = nx.get_node_attributes(G2, 'posxy')


    nx.draw_networkx_edges(G2, positions, node_size=0, edge_color='k',
                     font_size=0,  width=1)


    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    #limits=plt.axis('off')
    plt.draw()


def draw_edges_label(mesh, name, dpi):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]
    plt.figure(name, dpi=dpi)
    G = nx.Graph()


    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    bound_middle = {}
    iant = mesh.boundary_nodes[0, 0]

    cont = 0
    for i,e1,e2 in mesh.boundary_nodes:
        if i == iant:
            cont += 1
            bound_middle[i] = cont
        else:
            cont = 1
        iant = i

    edge_labels = {}
    cont = 0
    for i,e1,e2 in mesh.boundary_nodes:
        cont += 1
        if cont == int(bound_middle[i]/2.0):
            edge_labels[e1, e2] = str(i)
        if cont == bound_middle[i]:
            cont = 0

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edge_labels(G, positions, edge_labels, label_pos=0.5,
                                 font_size=9, font_color='b')

    plt.axes().set_aspect('equal')
    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    plt.draw()



def draw_nodes_label(mesh, name, dpi):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]
    plt.figure(name, dpi=dpi)
    G = nx.Graph()

    label = {}
    for i in range(len(X)):
        label[i] = i
        G.add_node(i, posxy=(X[i], Y[i]))


    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_nodes(G, positions, node_color='w', node_size=140,
                           node_shape='s')
    nx.draw_networkx_labels(G,positions,label,font_size=9)
    plt.axes().set_aspect('equal')
    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    plt.draw()



def draw_domain(mesh, name, dpi, color):

    c = mesh.nodes_coord

    plt.figure(name, dpi=dpi, frameon=False)

    X, Y = c[:, 0], c[:, 1]

    G = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G.add_node(i, posxy=(X[i], Y[i]))

    for plTag, lineTag in mesh.physicalLine.items():
        lineNodes = mesh.line[lineTag]

        G.add_edge(lineNodes[0], lineNodes[1])

    positions = nx.get_node_attributes(G, 'posxy')

    nx.draw_networkx_edges(G, positions, edge_color=color,
                     font_size=0, width=1, origin='lower')

    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)

    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    # limits=plt.axis('off')
    plt.draw()



def draw_elements(mesh, name, dpi, color):

    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]
    plt.figure(name, dpi=dpi, frameon = False)
    G2 = nx.Graph()

    label = []
    for i in range(len(X)):
        label.append(i)
        G2.add_node(i, posxy=(X[i], Y[i]))

    if mesh.gmsh == 1.0:
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1


    for i in range(len(mesh.ele_conn)):
        G2.add_cycle([mesh.ele_conn[i, 0],
                     mesh.ele_conn[i, 1],
                     mesh.ele_conn[i, 3],
                     mesh.ele_conn[i, 2]], )


    edge_line_nodes = {}
    for i in range(len(mesh.boundary_nodes[:, 0])):
        edge_line_nodes[(mesh.boundary_nodes[i, 1], mesh.boundary_nodes[i,
                                                                        2])] \
            = mesh.boundary_nodes[i, 0]


    positions = nx.get_node_attributes(G2, 'posxy')

    nx.draw_networkx(G2, positions, node_size=0, edge_color=color,
                     font_size=0,  width=1)


    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    #limits=plt.axis('off')


    plt.draw()

def draw_deformed_elements2(mesh, a, name, dpi, magf, color, G2):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]
    dX, dY = X + a[::2]*magf, Y + a[1::2]*magf


    if mesh.gmsh == 1.0:
        # correct nodes ordering
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1.0

    label2 = []
    for i in range(len(dX)):
        label2.append(i)
        G2.add_node(i, posxy2=(dX[i], dY[i]))

    for i in range(len(mesh.ele_conn)):
        G2.add_cycle([mesh.ele_conn[i, 0],
                     mesh.ele_conn[i, 1],
                     mesh.ele_conn[i, 3],
                     mesh.ele_conn[i, 2]], )

    positions2 = nx.get_node_attributes(G2, 'posxy2')

    frame = nx.draw_networkx_edges(G2, positions2, node_size=0, edge_color=color,
                     font_size=0, width=1, alpha=1)

    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    # limits=plt.axis('off')
    return frame

def draw_deformed_elements(mesh, a, name, dpi, magf, color, la):
    c = mesh.nodes_coord

    X, Y = c[:, 0], c[:, 1]
    dX, dY = X + a[::2]*magf, Y + a[1::2]*magf
    plt.figure(name, dpi=dpi, frameon=False)

    G2 = nx.Graph()

    if mesh.gmsh == 1.0:
        # correct nodes ordering
        temp = np.copy(mesh.ele_conn[:, 2])
        mesh.ele_conn[:, 2] = mesh.ele_conn[:, 3]
        mesh.ele_conn[:, 3] = temp
        mesh.gmsh += 1.0

    label2 = []
    for i in range(len(dX)):
        label2.append(i)
        G2.add_node(i, posxy2=(dX[i], dY[i]))

    for i in range(len(mesh.ele_conn)):
        G2.add_cycle([mesh.ele_conn[i, 0],
                     mesh.ele_conn[i, 1],
                     mesh.ele_conn[i, 3],
                     mesh.ele_conn[i, 2]], )

    positions2 = nx.get_node_attributes(G2, 'posxy2')

    nx.draw_networkx(G2, positions2, node_size=0, edge_color=color,
                     font_size=0, width=1, alpha=la)

    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    # limits=plt.axis('off')
    plt.draw()


def draw_deformed_domain(mesh, a, name, dpi, magf, color):
    c = mesh.nodes_coord

    plt.figure(name, dpi=dpi, frameon = False)
    bn = mesh.boundary_nodes

    adX = a[::2]
    adY = a[1::2]

    X, Y = c[bn[:, 1], 0], c[bn[:, 1], 1]
    dX, dY = c[bn[:, 1], 0] + adX[bn[:, 1]]*magf, c[bn[:, 1], 1] +  adY[bn[:,
                                                                        1]]*magf

    G2 = nx.Graph()

    label2 = []
    for i in range(len(dX)):
        label2.append(i)
        G2.add_node(i, posxy2=(dX[i], dY[i]))

    for i in range(len(bn[:, 0]) - 1):
       G2.add_edge(i, i+1)

    G2.add_edge(len(bn[:, 0]) - 1, 0)

    positions2 = nx.get_node_attributes(G2, 'posxy2')

    nx.draw_networkx(G2, positions2, node_size=0, edge_color=color,
                        font_size=0, width=1)

    plt.axes().set_aspect('equal')

    plt.axes().autoscale_view(True, True, True)
    plt.margins(y=0.1, x=0.1, tight=False)
    #limits=plt.axis('off')

    plt.draw()


def draw_deformed_domain2(mesh, a, name, dpi, magf, color, G2):
    c = mesh.nodes_coord

    bn = mesh.boundary_nodes

    adX = a[::2]
    adY = a[1::2]

    X, Y = c[bn[:, 1], 0], c[bn[:, 1], 1]
    dX, dY = c[bn[:, 1], 0] + adX[bn[:, 1]]*magf, c[bn[:, 1], 1] +  adY[bn[:,
                                                                        1]]*magf


    label2 = []
    for i in range(len(dX)):
        label2.append(i)
        G2.add_node(i, posxy2=(dX[i], dY[i]))

    for i in range(len(bn[:, 0]) - 1):
       G2.add_edge(i, i+1)

    G2.add_edge(len(bn[:, 0]) - 1, 0)

    positions2 = nx.get_node_attributes(G2, 'posxy2')

    frame = nx.draw_networkx_edges(G2, positions2, node_size=3, edge_color=color,
                        font_size=0, width=1)


    #plt.axes().set_aspect('equal')

    #plt.axes().autoscale_view(True, True, True)
    #plt.margins(y=0.1, x=0.1, tight=False)
    #limits=plt.axis('off')
    return frame





def tricontourf_deformed(a, mesh, d,  name, cmap, dpi, magf):
    """Plot contour with the tricoutour function and the boundary line with
    the boundary node.

    """
    fig = plt.figure(name, dpi=dpi)
    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    c = mesh.nodes_coord
    bn = mesh.boundary_nodes

    xx, yy, zz = c[:, 0]+d[::2]*magf, c[:, 1]+d[1::2]*magf, a

    dx = d[::2]*magf
    dy = d[1::2]*magf

    ccx = np.append(c[bn[:, 1], 0]+dx[bn[:, 1]], c[bn[0, 1], 0]+dx[bn[0, 1]])
    ccy = np.append(c[bn[:, 1], 1]+dy[bn[:, 1]], c[bn[0, 1], 1]+dy[bn[0, 1]])
    plt.plot(ccx, ccy, '-k')

    triangles = []
    for n1, n2, n3, n4 in mesh.ele_conn:
        triangles.append([n1, n2, n3])
        triangles.append([n1, n3, n4])

    triangles = np.asarray(triangles)

    CS2 = plt.tricontourf(xx, yy, triangles, zz, N=10, origin='lower',
                          cmap=cmap)

    CS3 = plt.tricontour(xx, yy, triangles, zz, N=10, origin='lower',colors='k')


    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size=0.3, pad=0.1)

    cbar = plt.colorbar(CS2, cax=cax)
    cbar.ax.set_label(name, fontsize=12)

    plt.clabel(CS3, fontsize=8, colors='k', fmt='%1.1f')

    limits=plt.axis('off')
    # plt.savefig('1.png', transparent=True, dpi=300)
    #plt.axes().set_aspect('equal')
    #plt.axes().autoscale_view(True, True, True)

    plt.draw()


def tricontourf(a, mesh, name, cmap, dpi):
    """Plot contour with the tricoutour function and the boundary line with
    the boundary node.

    """
    fig = plt.figure(name, dpi=dpi, frameon = False, facecolor='None')
    ax1 = fig.add_subplot(1, 1, 1, aspect='equal')
    c = mesh.nodes_coord
    bn = mesh.boundary_nodes

    xx, yy, zz = c[:, 0], c[:, 1], a

    ccx = np.append(c[bn[:, 1], 0], c[bn[0, 1], 0])
    ccy = np.append(c[bn[:, 1], 1], c[bn[0, 1], 1])

    triangles = []
    for n1, n2, n3, n4 in mesh.ele_conn:
        triangles.append([n1, n2, n3])
        triangles.append([n1, n3, n4])

    triangles = np.asarray(triangles)

    lev = 20

    try:
        CS2 = plt.tricontourf(xx, yy, triangles, zz, lev,
                              origin='lower',
                          cmap=cmap, antialiased=True)
    except ValueError:  #raised if `y` is empty.
        pass

    try:
        CS3 = plt.tricontour(xx, yy, triangles, zz, lev,
                             origin='lower',
                             colors='k')
        plt.clabel(CS3, fontsize=8, colors='k', fmt='%1.1f')
    except ValueError:
        pass

    plt.xticks(np.arange(min(xx), max(xx)+1, 1.0))
    plt.yticks(np.arange(min(yy), max(yy)+1, 1.0))

    #remove label and tick of axis
    plt.gca().axes.get_xaxis().set_visible(True)
    plt.gca().axes.get_yaxis().set_visible(True)

    #remove background and rectangular frame
    fig.patch.set_visible(False)
    ax1.patch.set_visible(False)
    ax1.axis('off')

    #plot a solid line in the boundary
    plt.plot(ccx , ccy, '-k')
    #plt.scatter(xx, yy, c=zz)
    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)

    #adjusts the contour color bar depending on the aspect ratio of domain
    if max(yy) < max(xx)/2.5:
        pos = 'bottom'
        pad = 0.6
    else:
        pos = 'right'
        pad = 0.3

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes(pos, size=0.2, pad=pad)

    if max(yy) < max(xx)/2.5:
        ori = 'horizontal'
    else:
        ori = 'vertical'

    cbar = plt.colorbar(CS2, cax=cax, orientation=ori)
    cbar.set_label(name, fontsize=12)


    # plt.savefig('1.png', transparent=True, dpi=300)
    #plt.axes().set_aspect('equal')
    #plt.axes().autoscale_view(True, True, True)

    plt.draw()


def tricontourf2(a, mesh, name, cmap, dpi):
    """Plot contour with the tricoutour function and the boundary line with
    the boundary node.

    """
    c = mesh.nodes_coord
    bn = mesh.boundary_nodes

    xx, yy, zz = c[:, 0], c[:, 1], a

    ccx = np.append(c[bn[:, 1], 0], c[bn[0, 1], 0])
    ccy = np.append(c[bn[:, 1], 1], c[bn[0, 1], 1])

    triangles = []
    for n1, n2, n3, n4 in mesh.ele_conn:
        triangles.append([n1, n2, n3])
        triangles.append([n1, n3, n4])

    triangles = np.asarray(triangles)

    lev = 20

    CS2 = plt.tricontourf(xx, yy, triangles, zz, lev,
                              origin='lower',
                          cmap=cmap, antialiased=True)

    '''

    #remove label and tick of axis
    plt.gca().axes.get_xaxis().set_visible(True)
    plt.gca().axes.get_yaxis().set_visible(True)

    #remove background and rectangular frame
    fig.patch.set_visible(False)
    ax1.patch.set_visible(False)
    ax1.axis('off')

    #plot a solid line in the boundary
    plt.plot(ccx , ccy, '-k')
    #plt.scatter(xx, yy, c=zz)
    plt.xlabel(r'$x$', fontsize=14)
    plt.ylabel(r'$y$', fontsize=14)

    #adjusts the contour color bar depending on the aspect ratio of domain
    if max(yy) < max(xx)/2.5:
        pos = 'bottom'
        pad = 0.6
    else:
        pos = 'right'
        pad = 0.3

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes(pos, size=0.2, pad=pad)

    if max(yy) < max(xx)/2.5:
        ori = 'horizontal'
    else:
        ori = 'vertical'
        '''
    #cbar = plt.colorbar(CS2, cax=cax, orientation=ori)
    #cbar.set_label(name, fontsize=12)


    # plt.savefig('1.png', transparent=True, dpi=300)
    #plt.axes().set_aspect('equal')
    #plt.axes().autoscale_view(True, True, True)

    #plt.draw()
    return CS2
