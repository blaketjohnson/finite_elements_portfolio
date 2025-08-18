import numpy as np

def jacobian_2D(xy_e, dN_dxi, dN_deta):
    """
    xy_e: (nen, 2) nodal (x,y) for one element
    dN_dxi, dN_deta: (nen,) parent-derivatives at (xi,eta)
    Returns J(2x2), detJ, invJ.
    """
    J = np.zeros((2,2))
    J[0,0] = (dN_dxi  @ xy_e[:,0])
    J[0,1] = (dN_dxi  @ xy_e[:,1])
    J[1,0] = (dN_deta @ xy_e[:,0])
    J[1,1] = (dN_deta @ xy_e[:,1])
    detJ = np.linalg.det(J)
    if detJ <= 0:
        raise ValueError("Non-positive detJ (check element mapping / node order)")
    invJ = np.linalg.inv(J)
    return J, detJ, invJ
