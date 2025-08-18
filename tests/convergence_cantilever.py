"""
Convergence study for cantilever tip displacement vs Euler-Bernoulli beam theory.
Saves docs/convergence_cantilever.png
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from fem.elements import K_structural_Q4
from fem.materials import D_plane_stress
from fem.assembly import assemble_global
from fem.bc import apply_dirichlet

def make_structured_Q4_mesh(nx, ny, Lx=1.0, Ly=0.2):
    xs = np.linspace(0, Lx, nx+1)
    ys = np.linspace(0, Ly, ny+1)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    xy = np.column_stack([X.ravel(), Y.ravel()])
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

def solve_disp(nx, ny, Ty=-1e5, t=0.01, Lx=1.0, Ly=0.2, E=210e9, nu=0.3):
    D = D_plane_stress(E, nu)
    xy, IEN, X, Y = make_structured_Q4_mesh(nx, ny, Lx=Lx, Ly=Ly)

    Ks = []
    for e in range(IEN.shape[0]):
        nodes = IEN[e]
        xy_e = xy[nodes]
        ke = K_structural_Q4(xy_e, D, t=t)
        Ks.append(ke)

    K = assemble_global(Ks, IEN, ndofs_per_node=2)
    f = np.zeros(K.shape[0])

    nxn = nx+1
    left_nodes = [j*(nxn) + 0 for j in range(ny+1)]
    clamp_dofs = np.array([n*2 for n in left_nodes] + [n*2+1 for n in left_nodes])

    right_nodes = [j*(nxn) + nx for j in range(ny+1)]
    for n in right_nodes:
        f[n*2+1] += Ty * (Ly / len(right_nodes)) * t

    apply_dirichlet(K, f, clamp_dofs, np.zeros_like(clamp_dofs, dtype=float))

    u = np.linalg.solve(K, f)
    U = u.reshape(-1,2)

    # Tip vertical displacement: mean Uy on right edge
    Uy_right = U[right_nodes, 1]
    return np.mean(Uy_right)

def main():
    # Beam theory tip deflection: delta = P L^3 / (3 E I)
    L, h, t = 1.0, 0.2, 0.01
    Ty = -1e5
    P = Ty * h * t
    E = 210e9
    I = t * h**3 / 12.0
    delta_ref = P * L**3 / (3 * E * I)
    print("Reference (beam) tip deflection:", delta_ref)

    ns = [4, 8, 12, 16, 24, 32]
    errs = []
    hs = []
    for nx in ns:
        ny = max(2, nx // 5)
        Uy_tip = solve_disp(nx, ny, Ty=Ty, t=t, Lx=L, Ly=h, E=E)
        err = abs(Uy_tip - (-delta_ref))
        errs.append(err)
        hs.append(L/nx)
        print(f"nx={nx:3d}, ny={ny:2d}, Uy_tip={Uy_tip:.6e}, err={err:.3e}")

    # Plot error vs element size h
    fig = plt.figure()
    plt.loglog(hs, errs, marker="o")
    plt.gca().invert_xaxis()
    plt.xlabel("Element size h (m)")
    plt.ylabel("Tip deflection error (m)")
    plt.title("Cantilever tip deflection convergence (Q4)")
    os.makedirs("docs", exist_ok=True)
    plt.savefig("docs/convergence_cantilever.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
