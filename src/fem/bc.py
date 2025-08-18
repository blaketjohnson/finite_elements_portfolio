import numpy as np

def apply_dirichlet(K, f, dof_indices, values):
    """
    Impose u[dof_indices] = values. Modifies K,f in-place and returns them.
    """
    for dof, val in zip(dof_indices, values):
        K[dof,:] = 0.0
        K[:,dof] = 0.0
        K[dof,dof] = 1.0
        f[dof] = val
    return K, f
