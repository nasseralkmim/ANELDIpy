__author__ = 'Nasser'

import numpy as np
import math


class Matrices:
    """Build the elemental matrices.

    Creates an object that has as attributes the elemental matrices.

    """
    def __init__(self, objectmesh):
        self.mesh = objectmesh

    def stiffness(self, material):
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
        mesh = self.mesh

        self.K = np.zeros((8, 8, mesh.num_ele))

        self.gp = mesh.chi / math.sqrt(3.0)

        for surf, mat in material.items():
            E = mat[0]
            nu = mat[1]

            C = np.zeros((3,3))

            C[0, 0] = 1.0
            C[1, 1] = 1.0
            C[1, 0] = nu
            C[0, 1] = nu
            C[2, 2] = (1.0 - nu)/2.0
            self.C = (E/(1.0-nu**2.0))*C

            for e in range(mesh.num_ele):
                if surf == mesh.ele_surface[e, 1]:
                    for w in range(4):
                        mesh.basisFunction2D(self.gp[w])
                        mesh.eleJacobian(mesh.nodes_coord[
                            mesh.ele_conn[e, :]])

                        if mesh.detJac <= 0.0:
                            print('Error: non-positive Jacobian - '
                                  'check element nodes numbering')
                            print('Element', e)

                        dp_xi = mesh.dphi_xi
                        B = np.array([
                            [dp_xi[0, 0], 0, dp_xi[0, 1], 0, dp_xi[0, 2], 0,
                             dp_xi[0, 3], 0]
                                     ,
                            [0, dp_xi[1, 0], 0, dp_xi[1, 1], 0, dp_xi[1, 2], 0,
                             dp_xi[1, 3]]
                                     ,
                            [dp_xi[1, 0], dp_xi[0, 0], dp_xi[1, 1], dp_xi[0, 1],
                            dp_xi[1, 2], dp_xi[0, 2],dp_xi[1, 3], dp_xi[0, 3]]])

                        CB = np.dot(self.C, B)

                        self.K[:, :, e] += (np.dot(np.transpose(B), CB) * mesh.detJac)



    def body_forces(self, q):
        """Build the load vector for the internal distributed load

        Args:
            mesh: object that includes the mesh attributes and methods for basis
                functions and transformation jacobians
            q (array of functions): internal distributed load

        """
        mesh = self.mesh
        self.P0q = np.zeros((8, mesh.num_ele))

        for e in range(mesh.num_ele):
            for w in range(4):
                mesh.basisFunction2D(self.gp[w])
                mesh.eleJacobian(mesh.nodes_coord[
                    mesh.ele_conn[e]])

                x1_o_e1e2, x2_o_e1e2 = mesh.mapping(e)

                load = q(x1_o_e1e2, x2_o_e1e2)

                self.P0q[0, e] += load[0]*mesh.phi[0]*mesh.detJac
                self.P0q[1, e] += load[1]*mesh.phi[0]*mesh.detJac
                self.P0q[2, e] += load[0]*mesh.phi[1]*mesh.detJac
                self.P0q[3, e] += load[1]*mesh.phi[1]*mesh.detJac
                self.P0q[4, e] += load[0]*mesh.phi[2]*mesh.detJac
                self.P0q[5, e] += load[1]*mesh.phi[2]*mesh.detJac
                self.P0q[6, e] += load[0]*mesh.phi[3]*mesh.detJac
                self.P0q[7, e] += load[1]*mesh.phi[3]*mesh.detJac

                
    def mass(self, material):
        """Build the elemental  mass matrix 8x8

        """
        mesh = self.mesh
        self.M = np.zeros((8, 8, mesh.num_ele))

        for surf, mat in material.items():
            rho = mat[2]

            for e in range(mesh.num_ele):
                if surf == mesh.ele_surface[e, 1]:
                    for w in range(4):
                        mesh.basisFunction2D(self.gp[w])
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

                        Jac = mesh.detJac
                        self.M[;, ;, e] += np.dot(np.transpose(N), N) * rho * Jac
        
        
                

   
