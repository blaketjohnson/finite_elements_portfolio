import numpy as np
from .shapes import shape_Q4
from .jacobian import jacobian_2D

def q4_B(xy_e, xi=0.0, eta=0.0):
    """Return B (3x8), detJ at (xi,eta) for a Q4 element."""
    N, dN_dxi, dN_deta = shape_Q4(xi, eta)
    J, detJ, invJ = jacobian_2D(xy_e, dN_dxi, dN_deta)
    dN_dx = invJ[0,0]*dN_dxi + invJ[0,1]*dN_deta
    dN_dy = invJ[1,0]*dN_dxi + invJ[1,1]*dN_deta
    B = np.zeros((3,8))
    for a in range(4):
        B[0, 2*a  ] = dN_dx[a]
        B[1, 2*a+1] = dN_dy[a]
        B[2, 2*a  ] = dN_dy[a]
        B[2, 2*a+1] = dN_dx[a]
    return B, detJ

def von_mises_plane_stress(sig):
    """sig = [sx, sy, txy] -> von Mises in plane stress."""
    sx, sy, txy = sig
    return np.sqrt(sx**2 - sx*sy + sy**2 + 3.0*txy**2)

def nodal_von_mises_Q4(xy, IEN, D, u):
    """Compute an approximate nodal von Mises by averaging element-center values."""
    Nn = xy.shape[0]
    vm = np.zeros(Nn)
    cnt = np.zeros(Nn)

    for e, nodes in enumerate(IEN):
        xy_e = xy[nodes]
        edofs = []
        for a in nodes:
            edofs.extend([2*a, 2*a+1])
        ue = u[edofs]
        B, _ = q4_B(xy_e, 0.0, 0.0)  # center point
        eps = B @ ue  # [ex, ey, gxy]
        sig = D @ eps
        vm_e = von_mises_plane_stress(sig)
        for a in nodes:
            vm[a] += vm_e
            cnt[a] += 1.0

    cnt[cnt == 0] = 1.0
    return vm / cnt
