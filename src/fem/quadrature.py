import numpy as np

def gauss_legendre_1D(n):
    if n == 1:
        return np.array([0.0]), np.array([2.0])
    if n == 2:
        a = 1/np.sqrt(3)
        return np.array([-a, a]), np.array([1.0, 1.0])
    if n == 3:
        a = np.sqrt(3/5)
        return np.array([-a, 0.0, a]), np.array([5/9, 8/9, 5/9])
    raise ValueError("Only n=1..3 supported")

def gauss_quad_2x2():
    """2x2 Gauss for Q4"""
    xi1d, w1d = gauss_legendre_1D(2)
    pts = []
    wts = []
    for i,x in enumerate(xi1d):
        for j,y in enumerate(xi1d):
            pts.append((x,y))
            wts.append(w1d[i]*w1d[j])
    return np.array(pts), np.array(wts)

def triangle_area_rule(order=1):
    """
    Simple area rules on reference triangle (xi>=0, eta>=0, xi+eta<=1).
    order=1: one-point rule at (1/3,1/3) with weight 1/2
    order=3: three-point rule at (1/6,1/6),(2/3,1/6),(1/6,2/3), each weight 1/6
    """
    if order == 1:
        return np.array([[1/3, 1/3]]), np.array([0.5])
    if order == 3:
        pts = np.array([[1/6, 1/6],[2/3, 1/6],[1/6, 2/3]])
        wts = np.array([1/6, 1/6, 1/6])
        return pts, wts
    raise ValueError("triangle_area_rule supports order=1 or 3")
