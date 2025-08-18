import numpy as np
from .shapes import shape_Q4, shape_T3
from .jacobian import jacobian_2D
from .quadrature import gauss_quad_2x2, triangle_area_rule

def K_structural_Q4(xy_e, D, t=1.0):
    """
    Plane stress/strain Q4, 2x2 Gauss.
    xy_e: (4,2), D: (3,3), thickness t.
    """
    K = np.zeros((8,8))
    pts, wts = gauss_quad_2x2()
    for (xi,eta), w in zip(pts, wts):
        N, dN_dxi, dN_deta = shape_Q4(xi, eta)
        J, detJ, invJ = jacobian_2D(xy_e, dN_dxi, dN_deta)
        dN_dx  = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
        dN_dy  = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
        B = np.zeros((3,8))
        for a in range(4):
            B[0, 2*a  ] = dN_dx[a]
            B[1, 2*a+1] = dN_dy[a]
            B[2, 2*a  ] = dN_dy[a]
            B[2, 2*a+1] = dN_dx[a]
        K += (B.T @ D @ B) * detJ * t * w
    return K

def K_structural_T3(xy_e, D, t=1.0, order=1):
    K = np.zeros((6,6))
    pts, wts = triangle_area_rule(order=order)
    for (xi,eta), w in zip(pts, wts):
        N, dN_dxi, dN_deta = shape_T3(xi, eta)
        J, detJ, invJ = jacobian_2D(xy_e, dN_dxi, dN_deta)
        dN_dx  = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
        dN_dy  = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
        B = np.zeros((3,6))
        for a in range(3):
            B[0, 2*a  ] = dN_dx[a]
            B[1, 2*a+1] = dN_dy[a]
            B[2, 2*a  ] = dN_dy[a]
            B[2, 2*a+1] = dN_dx[a]
        K += (B.T @ D @ B) * detJ * t * w
    return K

def K_conduction_Q4(xy_e, k):
    K = np.zeros((4,4))
    pts, wts = gauss_quad_2x2()
    for (xi,eta), w in zip(pts, wts):
        N, dN_dxi, dN_deta = shape_Q4(xi, eta)
        J, detJ, invJ = jacobian_2D(xy_e, dN_dxi, dN_deta)
        dN_dx  = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
        dN_dy  = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
        gradN = np.vstack((dN_dx, dN_dy))  # (2,4)
        K += (gradN.T @ (k* np.eye(2)) @ gradN) * detJ * w
    return K

def K_conduction_T3(xy_e, k, order=1):
    K = np.zeros((3,3))
    pts, wts = triangle_area_rule(order=order)
    for (xi,eta), w in zip(pts, wts):
        N, dN_dxi, dN_deta = shape_T3(xi, eta)
        J, detJ, invJ = jacobian_2D(xy_e, dN_dxi, dN_deta)
        dN_dx  = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
        dN_dy  = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
        gradN = np.vstack((dN_dx, dN_dy))
        K += (gradN.T @ (k* np.eye(2)) @ gradN) * detJ * w
    return K

def K_bar_2node(x1, x2, E, A):
    """
    1D 2-node bar element stiffness on [x1,x2].
    K_e = (E*A/L) [[1, -1], [-1, 1]]
    """
    L = float(x2 - x1)
    if L <= 0:
        raise ValueError("Bar length must be positive")
    k = (E*A)/L
    return k * np.array([[1.0, -1.0], [-1.0, 1.0]])

