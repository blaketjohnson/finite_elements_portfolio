import numpy as np

def assemble_global(Ks, IEN, ndofs_per_node):
    """
    Ks: list of element stiffness matrices (each (nen*ndofs, nen*ndofs))
    IEN: (nelems, nen) array of element -> node connectivity (global node ids 0..N-1)
    ndofs_per_node: 1 (conduction) or 2 (ux,uy) for structural
    Returns K_global (N*ndofs x N*ndofs)
    """
    nelems, nen = IEN.shape
    max_node = int(np.max(IEN))
    Nnodes = max_node + 1
    K = np.zeros((Nnodes*ndofs_per_node, Nnodes*ndofs_per_node))

    for e in range(nelems):
        ke = Ks[e]
        nodes = IEN[e]
        edofs = []
        for a in nodes:
            for d in range(ndofs_per_node):
                edofs.append(a*ndofs_per_node + d)
        edofs = np.array(edofs, dtype=int)
        K[np.ix_(edofs, edofs)] += ke
    return K

def assemble_force_RHS(Fs, IEN, ndofs_per_node):
    nelems, nen = IEN.shape
    max_node = int(np.max(IEN))
    Nnodes = max_node + 1
    f = np.zeros(Nnodes*ndofs_per_node)
    for e in range(nelems):
        fe = Fs[e]
        nodes = IEN[e]
        edofs = []
        for a in nodes:
            for d in range(ndofs_per_node):
                edofs.append(a*ndofs_per_node + d)
        edofs = np.array(edofs, dtype=int)
        f[edofs] += fe
    return f
