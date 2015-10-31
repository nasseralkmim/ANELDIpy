import numpy as np
import math


def stiffness(mesh, material):
    """Build the elemental stiffness matrix.

    Runs over each individual element properties when the object methods
    are called. The object is the mesh and its methods are basisFunction2D
    which defines the phi function and its derivative;  eleJacobian which
    produce a Jacobian matrix and its determinant and also the derivative
    of the basis function with respect the spatial coordinates.

    .. note::

        How it works:

        1. Loop over elements.

        2. loop over the 4 possible combination of 2 GP for each direction.

        3. Call the methods to create the derivative of shape functions  and Jacobian for each GP combination.

        4. Build the stiffness matrix from a matrix multiplication.

    Gauss points from natural nodal coordinates::

         gp = [ [-1, -1]
                [ 1, -1]
                [ 1,  1]
                [-1,  1] ]/ sqrt(3)


    Args:
        mesh: object that includes the mesh atributes and methods for basis
            functions and transformation jacobians.
        k (function) : material properties for the constitutive relation.

    Return:
        K_ele (array of float): 3nd order array 4x4 with the elemental
        stiffness matrix.

    """

    k_ele = np.zeros((8, 8, mesh.num_ele))

    gp = mesh.chi / math.sqrt(3.0)

    for surf, mat in material.items():
        E = mat[0]
        nu = mat[1]

        c1 = np.zeros((3, 3))

        c1[0, 0] = 1.0
        c1[1, 1] = 1.0
        c1[1, 0] = nu
        c1[0, 1] = nu
        c1[2, 2] = (1.0 - nu)/2.0
        c = (E/(1.0-nu**2.0))*c1

        for e in range(mesh.num_ele):
            if surf == mesh.ele_surface[e, 1]:
                for w in range(4):
                    mesh.basisFunction2D(gp[w])
                    mesh.eleJacobian(mesh.nodes_coord[
                        mesh.ele_conn[e, :]])

                    if mesh.detJac <= 0.0:
                        print('Error: non-positive Jacobian - '
                              'check element nodes numbering')
                        print('Element', e)

                    dp_xi = mesh.dphi_xi
                    b = np.array([
                        [dp_xi[0, 0], 0, dp_xi[0, 1], 0, dp_xi[0,2], 0,
                         dp_xi[0, 3], 0],
                        [0, dp_xi[1, 0], 0, dp_xi[1, 1], 0, dp_xi[1, 2], 0,
                         dp_xi[1, 3]],
                        [dp_xi[1, 0], dp_xi[0, 0], dp_xi[1, 1], dp_xi[0, 1],
                         dp_xi[1, 2], dp_xi[0, 2],dp_xi[1, 3], dp_xi[0, 3]]])

                    cb = np.dot(c, b)

                    k_ele[:, :, e] += (np.dot(np.transpose(b),
                                              cb) * mesh.detJac)


def body_forces(mesh, q):
    """Build the load vector for the internal distributed load

    Args:
        mesh: object that includes the mesh attributes and methods for basis
            functions and transformation jacobians
        q (array of functions): internal distributed load

    """
    p0q_ele = np.zeros((8, mesh.num_ele))

    gp = mesh.chi / math.sqrt(3.0)

    for e in range(mesh.num_ele):
        for w in range(4):
            mesh.basisFunction2D(gp[w])
            mesh.eleJacobian(mesh.nodes_coord[
                mesh.ele_conn[e]])

            x1_o_e1e2, x2_o_e1e2 = mesh.mapping(e)

            load = q(x1_o_e1e2, x2_o_e1e2)

            p0q_ele[0, e] += load[0]*mesh.phi[0]*mesh.detJac
            p0q_ele[1, e] += load[1]*mesh.phi[0]*mesh.detJac
            p0q_ele[2, e] += load[0]*mesh.phi[1]*mesh.detJac
            p0q_ele[3, e] += load[1]*mesh.phi[1]*mesh.detJac
            p0q_ele[4, e] += load[0]*mesh.phi[2]*mesh.detJac
            p0q_ele[5, e] += load[1]*mesh.phi[2]*mesh.detJac
            p0q_ele[6, e] += load[0]*mesh.phi[3]*mesh.detJac
            p0q_ele[7, e] += load[1]*mesh.phi[3]*mesh.detJac

    return p0q_ele


def mass(mesh, material):
    """Build the elemental  mass matrix 8x8

    """
    m_ele = np.zeros((8, 8, mesh.num_ele))

    gp = mesh.chi / math.sqrt(3.0)

    for surf, mat in material.items():
        rho = mat[2]

        for e in range(mesh.num_ele):
            if surf == mesh.ele_surface[e, 1]:
                for w in range(4):
                    mesh.basisFunction2D(gp[w])
                    mesh.eleJacobian(mesh.nodes_coord[
                        mesh.ele_conn[e, :]])

                    if mesh.detJac <= 0.0:
                        print('Error: non-positive Jacobian - '
                              'check element nodes numbering')
                        print('Element', e)

                    n = mesh.phi

                    N = np.array([
                        [n[0], 0., n[1], 0., n[2], 0., n[3], 0.],
                        [0., n[0], 0., n[1], 0., n[2], 0., n[3]]
                    ])

                    jac = mesh.detJac
                    m_ele[:, :, e] += np.dot(np.transpose(N), N) * rho *\
                                       jac
    return m_ele