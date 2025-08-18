import numpy as np

def shape_Q4(xi, eta):
    """
    Bilinear quad (Q4) shape functions and derivatives in parent coords.
    Returns N (4,), dN_dxi (4,), dN_deta (4,).
    """
    N1 = 0.25*(1 - xi)*(1 - eta)
    N2 = 0.25*(1 + xi)*(1 - eta)
    N3 = 0.25*(1 + xi)*(1 + eta)
    N4 = 0.25*(1 - xi)*(1 + eta)
    N = np.array([N1, N2, N3, N4])

    dN_dxi  = 0.25*np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    dN_deta = 0.25*np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return N, dN_dxi, dN_deta

def shape_T3(xi, eta):
    """
    Linear triangle (T3) on area coordinates: xi >=0, eta >=0, xi+eta<=1.
    N = [1 - xi - eta, xi, eta]
    dN_dxi = [-1, 1, 0], dN_deta = [-1, 0, 1]
    """
    N = np.array([1 - xi - eta, xi, eta])
    dN_dxi  = np.array([-1.0, 1.0, 0.0])
    dN_deta = np.array([-1.0, 0.0, 1.0])
    return N, dN_dxi, dN_deta

def shape_T6(xi, eta):
    """
    Quadratic triangle (T6) used in your HW 9.1.
    Node order: (1:(0,0)), (2:(1,0)), (3:(0,1)), mids (4,5,6)
    Returns N(6,), dN_dxi(6,), dN_deta(6,).
    """
    L1 = 1 - xi - eta
    L2 = xi
    L3 = eta
    N1 = L1*(2*L1 - 1)
    N2 = L2*(2*L2 - 1)
    N3 = L3*(2*L3 - 1)
    N4 = 4*L2*L3
    N5 = 4*L3*L1
    N6 = 4*L1*L2
    N = np.array([N1, N2, N3, N4, N5, N6])

    dL1_dxi,  dL1_deta  = -1.0, -1.0
    dL2_dxi,  dL2_deta  =  1.0,  0.0
    dL3_dxi,  dL3_deta  =  0.0,  1.0

    dN1_dL1 = 4*L1 - 1
    dN2_dL2 = 4*L2 - 1
    dN3_dL3 = 4*L3 - 1
    dN4_dL2, dN4_dL3 = 4*L3, 4*L2
    dN5_dL3, dN5_dL1 = 4*L1, 4*L3
    dN6_dL1, dN6_dL2 = 4*L2, 4*L1

    dN_dxi = np.array([
        dN1_dL1*dL1_dxi,
        dN2_dL2*dL2_dxi,
        dN3_dL3*dL3_dxi,
        dN4_dL2*dL2_dxi + dN4_dL3*dL3_dxi,
        dN5_dL3*dL3_dxi + dN5_dL1*dL1_dxi,
        dN6_dL1*dL1_dxi + dN6_dL2*dL2_dxi
    ], dtype=float)

    dN_deta = np.array([
        dN1_dL1*dL1_deta,
        dN2_dL2*dL2_deta,
        dN3_dL3*dL3_deta,
        dN4_dL2*dL2_deta + dN4_dL3*dL3_deta,
        dN5_dL3*dL3_deta + dN5_dL1*dL1_deta,
        dN6_dL1*dL1_deta + dN6_dL2*dL2_deta
    ], dtype=float)

    return N, dN_dxi, dN_deta
