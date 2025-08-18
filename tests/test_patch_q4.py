import numpy as np
from fem.elements import K_structural_Q4
from fem.materials import D_plane_stress

def test_q4_patch_linear_recovery():
    xy_e = np.array([[0,0],[1,0],[1,1],[0,1]], dtype=float)
    D = D_plane_stress(210e3, 0.3)
    K = K_structural_Q4(xy_e, D, t=1.0)
    w, _ = np.linalg.eig((K+K.T)/2)
    assert (w >= -1e-10).all()
    assert np.isfinite(K).all()
