# Throwaway: per-base-point section dofs of the three fieldsplit_1 subspaces on a tiny
# extruded periodic slice, and set relations LineSmooth vs Vanka vs Star patches.
import sys
sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")
from firedrake import *
from firedrake.preconditioners.asm import (get_local_ises_indices, build_star_indices,
                                           build_vanka_indices)
import numpy as np
import utils

nx, nz, deg = 4, 3, 2
m = PeriodicIntervalMesh(nx, 1.0)
mesh = ExtrudedMesh(m, nz, layer_height=1.0/nz)
V2 = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
Z = V2 * Vy * Vb
dm = mesh.topology_dm
print("layers =", mesh.layers, " nz =", nz)
print("base DMPlex: cells", dm.getHeightStratum(0), "vertices", dm.getHeightStratum(1), "dim", dm.getDimension())
names = ["RT_e (u)", "DG1xDG1 (u_y)", "DG1xCG2 (b)"]
for name, W in zip(names, Z):
    sec = W.dm.getLocalSection()
    nodes_per_entity = mesh.make_dofs_per_plex_entity(W.finat_element.entity_dofs())
    print(f"{name}: nodes_per_entity (on, above) per base dim = {nodes_per_entity}, block_size={W.block_size}")
    for codim, lab in [(0, "cell"), (1, "vertex")]:
        ps = range(*dm.getHeightStratum(codim))
        print(f"   {lab}: section dofs per point = {sorted(set(sec.getDof(p) for p in ps))}")
sec_mixed = Z.dm.getDefaultSection()
print("mixed Z.dm.getDefaultSection() handle:", sec_mixed.handle, "(0 means NULL section)")

# LineSmooth-style patches, mixed-aware
loc = get_local_ises_indices(Z)
secs = [W.dm.getLocalSection() for W in Z]
def line_patch(p):
    idx = []
    for i, W in enumerate(Z):
        d = secs[i].getDof(p)
        if d <= 0:
            continue
        o = secs[i].getOffset(p)
        idx.extend(loc[i][o*W.block_size:(o+d)*W.block_size])
    idx = np.asarray(idx, dtype=PETSc.IntType)
    return set(idx[idx >= 0].tolist())
cells = list(range(*dm.getHeightStratum(0)))
verts = list(range(*dm.getHeightStratum(1)))
col = {c: line_patch(c) for c in cells}
fac = {v: line_patch(v) for v in verts}
print("cell-column patch sizes:", sorted(set(len(s) for s in col.values())), " formula 14nz+4 =", 14*nz+4)
print("facet-line patch sizes:", sorted(set(len(s) for s in fac.values())), " formula 2nz =", 2*nz)
allp = list(col.values()) + list(fac.values())
tot = sum(len(s) for s in allp)
union = set().union(*allp)
print("LineSmooth 0,1: sum of sizes", tot, " union", len(union), " total dofs", Z.dim(),
      " pairwise disjoint:", tot == len(union), " covers all:", len(union) == Z.dim())

# Firedrake's own Vanka (closure of base cell) and Star (star of base vertex)
Zsplit = (tuple(Z), ())
locsplit = (loc, ())
for c in cells:
    vk = set(build_vanka_indices(Zsplit, locsplit, dm, "natural", "", True, c).indices.tolist())
    cone = dm.getCone(c)
    expect = col[c] | fac[cone[0]] | fac[cone[1]]
    assert vk == expect, (c, len(vk), len(expect))
print("Vanka(cell) == cell-column  U  facet-line(v0) U facet-line(v1): True for all cells; size", len(vk), " formula 18nz+4 =", 18*nz+4)
for v in verts:
    st = set(build_star_indices(Z, loc, dm, "natural", "", v).indices.tolist())
    supp = dm.getSupport(v)
    expect = fac[v] | col[supp[0]] | col[supp[1]]
    assert st == expect, (v, len(st), len(expect))
print("Star(vertex) == facet-line U cell-column(left) U cell-column(right): True for all vertices; size", len(st), " formula 30nz+8 =", 30*nz+8)
