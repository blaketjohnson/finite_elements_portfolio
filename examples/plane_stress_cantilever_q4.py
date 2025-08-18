"""
Q4 plane‑stress cantilever with traction on the right face.
Saves: docs/cantilever_deformed.png and docs/cantilever_disp_contour.png and docs/cantilever_vm.png
"""
import numpy as np
import os
import matplotlib.pyplot as plt
from fem.elements import K_structural_Q4
from fem.materials import D_plane_stress
from fem.assembly import assemble_global
from fem.bc import apply_dirichlet
from fem.post import nodal_von_mises_Q4

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

def main():
    nx, ny = 20, 4
    t = 0.01
    E, nu = 210e9, 0.3
    D = D_plane_stress(E, nu)

    xy, IEN, X, Y = make_structured_Q4_mesh(nx, ny, Lx=1.0, Ly=0.2)

    Ks = []
    for e in range(IEN.shape[0]):
        nodes = IEN[e]
        xy_e = xy[nodes]
        ke = K_structural_Q4(xy_e, D, t=t)
        Ks.append(ke)

    K = assemble_global(Ks, IEN, ndofs_per_node=2)
    f = np.zeros(K.shape[0])

    # Clamp x=0
    nxn = nx+1
    left_nodes = [j*(nxn) + 0 for j in range(ny+1)]
    clamp_dofs = np.array([n*2 for n in left_nodes] + [n*2+1 for n in left_nodes])

    # Tip traction Ty (downward) on right edge -> lumped nodal forces
    Ty = -1e5
    right_nodes = [j*(nxn) + nx for j in range(ny+1)]
    for n in right_nodes:
        f[n*2+1] += Ty * (0.2 / len(right_nodes)) * t

    apply_dirichlet(K, f, clamp_dofs, np.zeros_like(clamp_dofs, dtype=float))

    u = np.linalg.solve(K, f)
    U = u.reshape(-1,2)
    disp_mag = np.linalg.norm(U, axis=1)
    print("Displacement range |u|:", disp_mag.min(), disp_mag.max())

    # ----- Save deformed mesh overlay -----
    os.makedirs("docs", exist_ok=True)
    scale = 1000
    xy_def = xy + scale * U
    fig = plt.figure()
    plt.scatter(xy[:,0], xy[:,1], s=6, label="original")
    plt.scatter(xy_def[:,0], xy_def[:,1], s=6, label=f"deformed x{scale}")
    for e in IEN:
        loop = e.tolist() + [e[0]]
        pts0 = xy[loop]
        ptsd = xy_def[loop]
        plt.plot(pts0[:,0], pts0[:,1], linewidth=0.5)
        plt.plot(ptsd[:,0], ptsd[:,1], linewidth=0.5)
    plt.axis("equal"); plt.legend(); plt.title("Cantilever: original vs deformed")
    plt.tight_layout()
    plt.savefig("docs/cantilever_deformed.png", dpi=300, bbox_inches="tight")

    # ----- Vertical displacement contour (structured grid) -----
    Uy = U[:,1].reshape((ny+1, nx+1))
    fig = plt.figure()
    plt.pcolormesh(X, Y, Uy, shading="gouraud")
    plt.colorbar(label="Vertical displacement (m)")
    plt.xlabel("x"); plt.ylabel("y"); plt.title("Cantilever Uy contour")
    plt.savefig("docs/cantilever_disp_contour.png", dpi=300, bbox_inches="tight")

    # ----- Von Mises (nodal, averaged from element centers) -----
    vm = nodal_von_mises_Q4(xy, IEN, D, u).reshape((ny+1, nx+1))
    fig = plt.figure()
    plt.pcolormesh(X, Y, vm, shading="gouraud")
    plt.colorbar(label="von Mises (Pa)")
    plt.xlabel("x"); plt.ylabel("y"); plt.title("Cantilever von Mises (avg. from centers)")
    plt.savefig("docs/cantilever_vm.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
