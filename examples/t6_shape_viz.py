"""
Visualize T6 quadratic triangle shape functions over the parent triangle.
Saves docs/t6_shape_functions.png
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from fem.shapes import shape_T6

def main():
    # Sample points in parent triangle (xi>=0, eta>=0, xi+eta<=1)
    n = 80
    xis = np.linspace(0, 1, n)
    etas = np.linspace(0, 1, n)
    P = []
    for xi in xis:
        for eta in etas:
            if xi >= 0 and eta >= 0 and xi + eta <= 1:
                P.append((xi, eta))
    P = np.array(P)
    # fake triangulation in parametric space using Delaunay
    tri = mtri.Triangulation(P[:,0], P[:,1])
    # Evaluate shape functions
    Nvals = []
    for xi, eta in P:
        N, _, _ = shape_T6(xi, eta)
        Nvals.append(N)
    Nvals = np.array(Nvals)  # (M,6)

    # Plot 6 subplots
    fig, axes = plt.subplots(2, 3, figsize=(9,6))
    axes = axes.ravel()
    for i in range(6):
        ax = axes[i]
        tcf = ax.tricontourf(tri, Nvals[:,i], levels=20)
        fig.colorbar(tcf, ax=ax)
        ax.set_title(f"N{i+1}(xi,eta)")
        ax.set_xlabel("xi"); ax.set_ylabel("eta")
        ax.set_aspect("equal")
    fig.suptitle("T6 shape functions")
    os.makedirs("docs", exist_ok=True)
    plt.tight_layout()
    plt.savefig("docs/t6_shape_functions.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
