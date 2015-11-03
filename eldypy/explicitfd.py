import numpy as np

def scheme(T, N, u0, v0, m, p0, k):
    """ Perform a step on the finite difference scheme

    """
    dt = T/N
    u1 = u0 + dt*v0

    # TODO: transform u0 and u1 as vectors
    
    for i in range(N):
        # set variables to compute next step
        u_n1 = u1
        u_n0 = u0
        # compute next step
        a_n = np.linalg.solve(m, p0 - np.dot(k, u_n1))
        u_n2 = dt**2.0*a_n + 2*u_n1 - u_n0
        # update initial variables
        u_n0 = u_n1
        u_n1 = u_n2

        
