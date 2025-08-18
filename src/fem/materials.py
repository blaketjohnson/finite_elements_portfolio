import numpy as np

def D_plane_stress(E, nu):
    c = E / (1 - nu**2)
    return c * np.array([[1,    nu,   0],
                         [nu,   1,    0],
                         [0,    0,  (1-nu)/2]])

def D_plane_strain(E, nu):
    c = E / ((1+nu)*(1-2*nu))
    return c * np.array([[1-nu,   nu,        0],
                         [nu,     1-nu,      0],
                         [0,      0,   0.5- nu]])

def k_isotropic(k):
    """ scalar conductivity -> used as k * I (but assembly uses grad N, so we just pass k) """
    return float(k)
