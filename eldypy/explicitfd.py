import numpy as np
from eldypy import plotter


def scheme(u0, u1, m, p0, k, dt):
    """ Perform a step on the finite difference scheme

    """
    # set variables to compute next step
    u_n1 = u1
    u_n0 = u0
    # compute next step
    a_n = np.linalg.solve(m, p0 - np.dot(k, u_n1))
    u_n2 = dt**2.0*a_n + 2*u_n1 - u_n0
    return u_n2
