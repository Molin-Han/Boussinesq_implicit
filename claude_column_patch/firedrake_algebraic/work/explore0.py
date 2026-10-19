# exploration only: mesh, spaces, patches
import sys
sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")
from firedrake import *
import utils
import numpy as np
from firedrake.preconditioners.asm import (build_star_indices, build_vanka_indices,
                                           get_local_ises_indices, get_colors, ASMStarPC, ASMVankaPC)
nx, nz, length, height, refinement, deg = 5, 20, 3.0e5, 4.0e3, 2, 2
distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (refinement+1), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]
print("cells", mesh.num_cells(), "layers", mesh.layers, "geometric dim", mesh.geometric_dimension)
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg-1)
Pressure = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure
W1 = V_2D * Vy * Vb
print("dims", [V.dim() for V in W])
dm = mesh.topology_dm
print("dm dim", dm.getDimension(), "chart", dm.getChart(), "depth strata", [dm.getDepthStratum(d) for d in range(dm.getDimension()+1)])
ises_idx = get_local_ises_indices(W1)
for kind in ["star", "vanka"]:
    if kind == "star":
        colors = get_colors(dm, False, 0, distance=1)
        ises = [build_star_indices(W1, ises_idx, dm, "natural", "x_", c) for c in colors]
    else:
        colors = get_colors(dm, False, 1, distance=3)
        Z = (tuple(W1[i] for i in range(3)), tuple())
        Zi = (tuple(ises_idx[i] for i in range(3)), tuple())
        ises = [build_vanka_indices(Z, Zi, dm, "natural", "x_", True, c) for c in colors]
    sizes = [i.getSize() for i in ises]
    print(kind, "npatches", len(ises), "sizes", set(sizes))
