"""
Visual patch test for a single Q4 element under a linear displacement field.
Saves docs/patch_q4_plot.png
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from fem.materials import D_plane_stress
from fem.elements import K_structural_Q4
from fem.assembly import assemble_global
from fem.bc import apply_dirichlet

def q4_unit_square_mesh():
    xy = np.array([[0,0],[1,0],[1,1],[0,1]], dtype=float)
    IEN = np.array([[0,1,2,3]], dtype=int)
    return xy, IEN

def main():
    xy, IEN = q4_unit_square_mesh()
    D = D_plane_stress(210e9, 0.3)
    t = 1.0

    # Assemble (single element is fine)
    ke = K_structural_Q4(xy[IEN[0]], D, t=t)
    K = assemble_global([ke], IEN, ndofs_per_node=2)
    f = np.zeros(8)

    # Exact linear displacement field: u=ax+b, v=cy+d
    a,b,c,d = 0.01, 0.0, -0.02, 0.0
    U_exact = np.zeros((4,2))
    for i,(x,y) in enumerate(xy):
        U_exact[i,0] = a*x + b
        U_exact[i,1] = c*y + d

    # Impose exact displacements on all nodes (Dirichlet) -> FEM should reproduce exactly
    dofs = []
    vals = []
    for n in range(4):
        dofs.extend([2*n, 2*n+1])
        vals.extend([U_exact[n,0], U_exact[n,1]])
    apply_dirichlet(K, f, np.array(dofs), np.array(vals, dtype=float))

    u = np.linalg.solve(K, f).reshape(-1,2)
    err = np.linalg.norm(u - U_exact, axis=1)

    # Plot nodal error
    fig = plt.figure()
    plt.stem(range(4), err)   # old: plt.stem(range(4), err, use_line_collection=True)
    plt.xlabel("Node index")
    plt.ylabel("Error magnitude (m)")
    plt.title("Q4 Patch Test: nodal displacement error (should be ~1e-15)")
    os.makedirs("docs", exist_ok=True)
    plt.savefig("docs/patch_q4_plot.png", dpi=300, bbox_inches="tight")
    print("Max nodal error:", err.max())

if __name__ == "__main__":
    main()
