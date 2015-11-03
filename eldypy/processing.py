import numpy as np


def principal_stress_max(s11, s22, s12):
    """

    :param s11:
    :param s22:
    :param s12:
    :return:
    """
    sp_max = np.zeros(len(s11))
    for i in range(len(s11)):
        sp_max[i] = (s11[i]+s22[i])/2.0 + np.sqrt((s11[i] - s22[i])**2.0/2.0 +
                                               s12[i]**2.0)

    return sp_max


def principal_stress_min(s11, s22, s12):
    """

    :param s11:
    :param s22:
    :param s12:
    :return:
    """
    sp_min = np.zeros(len(s11))
    for i in range(len(s11)):
        sp_min[i] = (s11[i]+s22[i])/2. - np.sqrt((s11[i] - s22[i])**2./2. +
                                               s12[i]**2.)

    return sp_min


def stress_recovery(mesh, u, material):
    """Recover stress and strain from nodal displacement value

    :param mesh:
    :param d:
    :param C:
    :return:
    """
    stress_ele = np.zeros((3, 4, mesh.num_ele))
    strain_ele = np.zeros((3, 4, mesh.num_ele))

    for surf, mat in material.items():
        C = linear_elastic_constitutive(mat)

        for e in range(mesh.num_ele):
            if surf == mesh.ele_surface[e, 1]:
                for w in range(4):
                    # Evaluate basis function at nodes natural coordinates chi
                    mesh.basisFunction2D(mesh.chi[w]/np.sqrt(3.0))
                    mesh.eleJacobian(mesh.nodes_coord[
                        mesh.ele_conn[e, :]])

                    B = np.array([
                        [mesh.dphi_xi[0, 0], 0, mesh.dphi_xi[0, 1], 0,
                        mesh.dphi_xi[0, 2], 0, mesh.dphi_xi[0, 3], 0]
                                 ,
                        [0, mesh.dphi_xi[1, 0], 0, mesh.dphi_xi[1, 1], 0,
                        mesh.dphi_xi[1, 2], 0, mesh.dphi_xi[1, 3]]
                                 ,
                        [mesh.dphi_xi[1, 0], mesh.dphi_xi[0, 0],
                        mesh.dphi_xi[1, 1], mesh.dphi_xi[0, 1],
                        mesh.dphi_xi[1, 2], mesh.dphi_xi[0, 2],
                        mesh.dphi_xi[1, 3], mesh.dphi_xi[0, 3]]])

                    u_ele = np.array([
                        u[2*mesh.ele_conn[e, 0]],
                        u[2*mesh.ele_conn[e, 0]+1],
                        u[2*mesh.ele_conn[e, 1]],
                        u[2*mesh.ele_conn[e, 1]+1],
                        u[2*mesh.ele_conn[e, 2]],
                        u[2*mesh.ele_conn[e, 2]+1],
                        u[2*mesh.ele_conn[e, 3]],
                        u[2*mesh.ele_conn[e, 3]+1],
                    ])

                    e_ele = np.dot(B, u_ele)

                    strain_ele[:, w, e] = e_ele

                    # w represents each node
                    stress_ele[:, w, e] = np.dot(C, e_ele)

    # Nodal stress Averaged by number of elements sharing a node
    s11node = global_average_vector(stress_ele[0, :, :], mesh)
    s22node = global_average_vector(stress_ele[1, :, :], mesh)
    s12node = global_average_vector(stress_ele[2, :, :], mesh)

    s11node = np.reshape(s11node, len(s11node))
    s22node = np.reshape(s22node, len(s22node))
    s12node = np.reshape(s12node, len(s12node))

    s_ele = np.array([element_stress(s11node, mesh),
                    element_stress(s22node, mesh),
                    element_stress(s12node, mesh)])

    e11node = global_average_vector(strain_ele[0, :, :], mesh)
    e22node = global_average_vector(strain_ele[1, :, :], mesh)
    e12node = global_average_vector(strain_ele[2, :, :], mesh)

    e11node = np.reshape(e11node, len(s11node))
    e22node = np.reshape(e22node, len(s22node))
    e12node = np.reshape(e12node, len(s12node))

    e_ele = np.array([element_stress(e11node, mesh),
                     element_stress(e22node, mesh),
                     element_stress(e12node, mesh)])

    return [s11node, s22node, s12node], s_ele, e_ele


def global_average_vector(v_ele, mesh):


    v = np.zeros((mesh.num_nodes, 1))

    for e in range(mesh.num_ele):
        for i in range(len(mesh.ele_conn[e])):
            v[mesh.ele_conn[e, i]] += v_ele[i, e]/((mesh.ele_conn ==
                                                    mesh.ele_conn[e, i]).sum())

    return v


def element_stress(s, mesh):
    """

    :param s:
    :param mesh:
    :return:
    """
    s_ele = np.zeros(mesh.num_ele)

    for e in range(mesh.num_ele):
        s_e = np.array([
                    s[mesh.ele_conn[e, 0]],
                    s[mesh.ele_conn[e, 1]],
                    s[mesh.ele_conn[e, 2]],
                    s[mesh.ele_conn[e, 3]],
                ])

        s_ele[e] = (s_e[0] + s_e[1] + s_e[2] + s_e[3])/4.0

    return s_ele



def sub_matrices(K, P, displacement):

    ie = []
    for l in displacement(1, 1).keys():
        if l[0] == 'nodes' or l[0] == 'node':
            for n in l[1:]:
                rx = displacement(1,1)[l][0]
                ry = displacement(1,1)[l][1]


                if rx != 'free':
                   ie.append([n, 0, rx])

                if ry != 'free':
                    ie.append([n, 1, ry])

    ide = []
    de = np.zeros(len(ie))
    for i in range(len(de)):
        de[i] = ie[i][2]

        if ie[i][1] == 0:
            ide.append(2*ie[i][0])
        if ie[i][1] == 1:
            ide.append(2*ie[i][0]+1)

    idf = []
    for i in range(len(P)):
        if i in ide:
            continue
        else:
            idf.append(i)


    Ke = np.zeros((len(ie), len(ie)))
    Ke = K[np.ix_(ide, ide)]
    Kf = K[np.ix_(idf, idf)]
    Kef = K[np.ix_(ide, idf)]
    Kfe = K[np.ix_(idf, ide)]
    Pf = P[np.ix_(idf)]
    Pe = P[np.ix_(ide)]

    return Ke, Kf, Kef, Kfe, de, Pf, Pe, ide, idf


def join_dof_values(de, df, ide, idf):

    U = np.zeros(len(ide)+len(idf))
    for i in range(len(ide)):
        U[ide[i]] = de[i]

    for i in range(len(idf)):
        U[idf[i]] = df[i]

    return U


def linear_elastic_constitutive(mat):
    E = mat[0]
    nu = mat[1]
    C = np.zeros((3,3))
    C[0, 0] = 1.0
    C[1, 1] = 1.0
    C[1, 0] = nu
    C[0, 1] = nu
    C[2, 2] = (1.0 - nu)/2.0
    C = (E/(1.0-nu**2.0))*C
    return C


