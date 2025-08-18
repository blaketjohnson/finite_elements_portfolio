"""
2D Heat conduction on a unit square using Q4 elements (homework‑style).
Dirichlet: T=10 on left & bottom, insulated on right & top.
Saves a contour image to docs/heat_q4_temperature.png
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from fem.elements import K_conduction_Q4
from fem.assembly import assemble_global
from fem.bc import apply_dirichlet

def make_structured_Q4_mesh(nx, ny, Lx=1.0, Ly=1.0):
    xs = np.linspace(0, Lx, nx+1)
    ys = np.linspace(0, Ly, ny+1)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    xy = np.column_stack([X.ravel(), Y.ravel()])  # (N,2)
    IEN = []
    def nid(i,j): return j*(nx+1) + i
    for j in range(ny):
        for i in range(nx):
            n1 = nid(i, j)
            n2 = nid(i+1, j)
            n3 = nid(i+1, j+1)
            n4 = nid(i, j+1)
            IEN.append([n1,n2,n3,n4])
    return xy, np.array(IEN, dtype=int), X, Y

def main():
    nx, ny = 20, 20
    k = 4.0  # conductivity
    xy, IEN, X, Y = make_structured_Q4_mesh(nx, ny)
    ne = IEN.shape[0]

    Ks = []
    for e in range(ne):
        nodes = IEN[e]
        xy_e = xy[nodes]
        ke = K_conduction_Q4(xy_e, k)
        Ks.append(ke)

    K = assemble_global(Ks, IEN, ndofs_per_node=1)
    f = np.zeros(K.shape[0])

    # Dirichlet: T=10 on x=0 and y=0 (left & bottom)
    nxn = nx+1
    left = [i*(nx+1) for i in range(ny+1)]
    bottom = list(range(nxn))
    bc_nodes = sorted(set(left + bottom))
    bc_dofs = np.array(bc_nodes)  # 1 dof per node
    bc_vals = np.full_like(bc_dofs, 10.0, dtype=float)

    apply_dirichlet(K, f, bc_dofs, bc_vals)

    T = np.linalg.solve(K, f)

    print("Solved temperature field, min/max:", T.min(), T.max())

    # ----- Save temperature contour -----
    Tgrid = T.reshape((ny+1, nx+1))
    fig = plt.figure()
    plt.pcolormesh(X, Y, Tgrid, shading="gouraud")
    plt.colorbar(label="Temperature")
    plt.xlabel("x"); plt.ylabel("y"); plt.title("Q4 heat plate (Dirichlet left/bottom)")
    os.makedirs("docs", exist_ok=True)
    plt.savefig("docs/heat_q4_temperature.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
