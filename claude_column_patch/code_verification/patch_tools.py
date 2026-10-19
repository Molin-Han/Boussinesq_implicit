"""
Shared tools of the verification scripts (02_patch_dofs.py, 03_operator_serial_vs_parallel.py).

Every DoF of the mixed space V = V_2D * Vy * Vb (the space of field 1 of the fieldsplit on one
multigrid level) is
  * classified: field (0 = in-plane velocity u, 1 = u_y, 2 = b) and entity kind
    (0 = vertical facet, 1 = horizontal facet, 2 = cell interior),
  * located: (x, z) of its node,
  * given a partition-independent key (field, 'c' or 'v', x of the base cell midpoint or of the
    base vertex, position of the DoF inside the column of that base entity).
All names that the solver script could rebind (it does `from firedrake import *`) are avoided:
this module has its own namespace.
"""
import numpy as np
import firedrake as fd

KIND = {0: 'vertical facet', 1: 'horizontal facet', 2: 'cell interior'}
FIELD = {0: 'u (in-plane, RT)', 1: 'u_y (DG)', 2: 'b (DG x CG)'}


def canon(x, L, nd=3):
    """x modulo the period L, rounded; the periodic seam is always reported as 0."""
    y = round(float(x) % L, nd)
    if abs(y - L) < 10.0 ** (-nd + 1):
        y = 0.0
    return y


def scalar_companions(V):
    """Scalar nodal spaces whose cell-local nodes are, in this order, the cell-local nodes of V."""
    mesh = V.mesh()
    el = V.ufl_element()
    if type(el).__name__ == 'EnrichedElement':     # RT_horiz + RT_vert: strip HDiv, keep the tensor products
        parts = []
        for sub in el._elements:
            assert type(sub).__name__ == 'HDivElement', sub
            parts.append(fd.FunctionSpace(mesh, sub._element))
        # the cell-local nodes of the enriched element must be those of the parts, one part after the
        # other: checked on the entity dofs
        ed = V.finat_element.entity_dofs()
        shift = 0
        merged = {}
        for S in parts:
            for dim, ents in S.finat_element.entity_dofs().items():
                for e, dofs in ents.items():
                    merged.setdefault(dim, {}).setdefault(e, []).extend(int(d) + shift for d in dofs)
            shift += S.cell_node_map().arity
        for dim, ents in ed.items():
            for e, dofs in ents.items():
                assert sorted(int(d) for d in dofs) == sorted(merged[dim][e]), (dim, e, dofs, merged[dim][e])
        return parts
    return [V]


def node_table(V):
    """For every local node (owned + halo) of V: kind (0/1/2), x, z, base point, position in the column."""
    mesh = V.mesh()
    dm = mesh.topology_dm
    nlayers = mesh.layers - 1
    cmap = V.cell_node_map()
    vals = cmap.values_with_halo
    off = cmap.offset
    ncell, nloc = vals.shape
    nnode = V.node_set.total_size

    ed = V.finat_element.entity_dofs()
    loc_kind = np.full(nloc, -1, dtype=int)
    for (dim, ents) in ed.items():
        for e, dofs in ents.items():
            for d in dofs:
                if dim == (0, 1):
                    loc_kind[d] = 0          # base vertex x vertical interval = vertical facet
                elif dim == (1, 0):
                    loc_kind[d] = 1          # base cell x vertical vertex = horizontal facet
                elif dim == (1, 1):
                    loc_kind[d] = 2
                else:
                    raise RuntimeError(f"DoF on a 2D vertex {dim}: not expected for these spaces")
    assert (loc_kind >= 0).all()

    comps = scalar_companions(V)
    assert sum(S.cell_node_map().arity for S in comps) == nloc
    tabs = []
    X = fd.SpatialCoordinate(mesh)
    for S in comps:
        fx = fd.Function(S).interpolate(X[0])
        fz = fd.Function(S).interpolate(X[2])
        tabs.append((fx.dat.data_ro_with_halos.copy(), fz.dat.data_ro_with_halos.copy(),
                     S.cell_node_map().values_with_halo, S.cell_node_map().offset))

    kind = np.full(nnode, -1, dtype=int)
    x = np.full(nnode, np.nan)
    z = np.full(nnode, np.nan)
    for c in range(ncell):
        for layer in range(nlayers):
            nodes = vals[c] + layer * off
            kind[nodes] = loc_kind
            j0 = 0
            for (fx, fz, sv, so) in tabs:
                sn = sv[c] + layer * so
                x[nodes[j0:j0 + len(sn)]] = fx[sn]
                z[nodes[j0:j0 + len(sn)]] = fz[sn]
                j0 += len(sn)
    assert (kind >= 0).all(), "some node is not reached from any local cell"

    sec = V.dm.getLocalSection()
    ps, pe = dm.getChart()
    base = np.full(nnode, -1, dtype=int)
    pos = np.full(nnode, -1, dtype=int)
    for p in range(ps, pe):
        ndof = sec.getDof(p)
        if ndof > 0:
            o = sec.getOffset(p)
            base[o:o + ndof] = p
            pos[o:o + ndof] = np.arange(ndof)
    assert (base >= 0).all()
    return kind, x, z, base, pos


def base_point_x(mesh, L):
    """x of every point of the base plex (cell midpoint / vertex), partition independent (mod L)."""
    dm = mesh.topology_dm
    ps, pe = dm.getChart()
    cs, ce = dm.getHeightStratum(0)
    vs, ve = dm.getDepthStratum(0)
    cf = mesh.coordinates
    cv = cf.function_space().cell_node_map().values_with_halo
    data = cf.dat.data_ro_with_halos
    cell_numbering = mesh._cell_numbering
    px = np.full(pe - ps, np.nan)
    iscell = np.zeros(pe - ps, dtype=bool)
    x_lo, x_hi = {}, {}
    for c in range(cs, ce):
        fc = cell_numbering.getOffset(c)
        xx = data[cv[fc], 0]
        lo, hi = xx.min(), xx.max()
        px[c - ps] = 0.5 * (lo + hi)
        iscell[c - ps] = True
        x_lo[c], x_hi[c] = lo, hi
    ends = {}
    for c in range(cs, ce):
        for v in dm.getCone(c):
            ends.setdefault(v, []).append(c)
    for v in range(vs, ve):
        cells = ends.get(v, [])
        if len(cells) == 2:
            a = {canon(x_lo[cells[0]], L), canon(x_hi[cells[0]], L)}
            b = {canon(x_lo[cells[1]], L), canon(x_hi[cells[1]], L)}
            common = a & b
            if len(common) == 1:
                px[v - ps] = common.pop()
                continue
        px[v - ps] = np.nan      # vertex seen from one local cell only (edge of the halo): resolved later
    return px, iscell


def describe_space(V, L):
    """Tables for the mixed space V handed to the smoother (all arrays in the local mixed numbering)."""
    mesh = V.mesh().unique() if hasattr(V.mesh(), 'unique') else V.mesh()
    px, iscell = base_point_x(mesh, L)
    dm = mesh.topology_dm
    ps, pe = dm.getChart()
    ghost = np.array([dm.getLabelValue("pyop2_ghost", p) != -1 for p in range(ps, pe)])
    local_ises = [iset.indices for iset in V.dof_dset.local_ises]
    ntot = sum(len(i) for i in local_ises)
    field = np.full(ntot, -1, dtype=int)
    kind = np.full(ntot, -1, dtype=int)
    x = np.full(ntot, np.nan)
    z = np.full(ntot, np.nan)
    base = np.full(ntot, -1, dtype=int)
    pos = np.full(ntot, -1, dtype=int)
    owned = np.zeros(ntot, dtype=bool)
    for i, W in enumerate(V):
        assert W.block_size == 1
        k, xx, zz, b, p = node_table(W)
        idx = local_ises[i]
        field[idx] = i
        kind[idx] = k
        x[idx] = xx
        z[idx] = zz
        base[idx] = b
        pos[idx] = p
        owned[idx[:W.dof_dset.size]] = True
    assert (field >= 0).all()
    # x of the base vertices that are seen from one local cell only: taken from their facet DoFs
    for p in range(ps, pe):
        if np.isnan(px[p - ps]):
            sel = (base == p) & (field == 0)
            if sel.any():
                px[p - ps] = sorted(set(canon(v, L) for v in x[sel]))[0]
    return dict(mesh=mesh, px=px, iscell=iscell, ghost=ghost, field=field, kind=kind, x=x, z=z, base=base,
                pos=pos, ntot=ntot, owned=owned, L=L)


def keys_of(T, idx):
    """Partition-independent keys of the DoFs idx (local mixed numbering)."""
    out = []
    for i in idx:
        p = T['base'][i]
        out.append((int(T['field'][i]), 'c' if T['iscell'][p] else 'v', canon(T['px'][p], T['L']), int(T['pos'][i])))
    return out


def key_arrays(T, idx):
    """The keys as arrays (field, iscell, bx, pos), for sorting and saving."""
    p = T['base'][idx]
    bx = np.array([canon(v, T['L']) for v in T['px'][p]])
    return T['field'][idx].astype(np.int64), T['iscell'][p].astype(np.int64), bx, T['pos'][idx].astype(np.int64)


def key_values(T, idx, which):
    """Deterministic test vectors defined by the keys only (so identical for every partition).

    which = 'smooth'  : smooth in x and z
            'rough'   : pseudo-random, O(1), different for every DoF
            'mixed'   : smooth + 0.3 * rough
    """
    f, c, bx, pos = key_arrays(T, idx)
    xx = T['x'][idx] % T['L'] / T['L']
    zz = T['z'][idx]
    zz = zz / max(zz.max(), 1.0) if len(zz) else zz
    smooth = np.sin(2 * np.pi * (bx / T['L']) * 2 + 0.3 * f) * np.cos(np.pi * pos / (pos.max() + 1.0) * 1.5 + 0.2 * c)
    h = np.sin(12.9898 * (bx / T['L'] * 977.0) + 78.233 * pos + 37.719 * f + 11.13 * c) * 43758.5453
    rough = 2.0 * (h - np.floor(h)) - 1.0
    if which == 'smooth':
        return smooth
    if which == 'rough':
        return rough
    return smooth + 0.3 * rough
