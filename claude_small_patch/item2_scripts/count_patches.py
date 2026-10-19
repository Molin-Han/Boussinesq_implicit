# Throwaway: enumerate candidate smoother patches with Firedrake's own get_patches and count
# divergence-free fields inside each patch.  Tiny serial mesh, no solver run.
import sys
sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")
import numpy as np
from firedrake import *
from firedrake.petsc import PETSc
from firedrake.preconditioners import asm
import utils

nx, nz, L, H, deg = 4, 3, 4.0, 1.0, 2
base = PeriodicIntervalMesh(nx, L)
mesh = ExtrudedMesh(base, nz, layer_height=H/nz)
V = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
Wm = V * Vy * Vb
print("dofs: u", V.dim(), "uy", Vy.dim(), "b", Vb.dim(), "mixed", Wm.dim())
print("expected per column: u 6Nz+2 =", 6*nz+2, " uy 4Nz =", 4*nz, " b 4Nz+2 =", 4*nz+2)

opts = PETSc.Options()
def patches(cls, space, prefix, **kw):
    for k_, v_ in kw.items():
        opts[prefix + k_] = v_
    obj = cls.__new__(cls)
    obj.prefix = prefix
    ises = obj.get_patches(space)
    for k_ in kw:
        opts.delValue(prefix + k_)
    return [np.array(i.indices) for i in ises]

def report(name, ps, expect=None):
    sizes = sorted(set(len(p) for p in ps))
    print(f"{name}: {len(ps)} patches, sizes {sizes}" + (f"  expected {expect}" if expect else ""))
    return ps

report("ASMStarPC dim0 (mixed)   [committed star]", patches(asm.ASMStarPC, Wm, "s0_", construct_dim=0), f"{30*nz+8}")
report("ASMStarPC dim1 (mixed)   [cell column, no lateral facets]", patches(asm.ASMStarPC, Wm, "s1_", construct_dim=1), f"{14*nz+4}")
report("ASMVankaPC dim1 (mixed)  [single column of the study]", patches(asm.ASMVankaPC, Wm, "v1_", construct_dim=1), f"{18*nz+4}")
report("ASMVankaPC dim0 exclude all, entity (mixed) [facet line]", patches(asm.ASMVankaPC, Wm, "v0_", construct_dim=0, exclude_subspaces="0,1,2", include_type="entity"), f"{2*nz}")
report("ASMExtrudedStarPC dim0 (mixed) [2D vertex star]", patches(asm.ASMExtrudedStarPC, Wm, "e0_", construct_dim=0), "52 interior / 28 top,bottom; nx*(nz+1) patches")
report("ASMExtrudedStarPC dim0 coloring (mixed)", patches(asm.ASMExtrudedStarPC, Wm, "e0c_", construct_dim=0, use_coloring=True))
try:
    patches(asm.ASMLinesmoothPC, Wm, "l_", codims="0,1")
except NotImplementedError as e:
    print("ASMLinesmoothPC on mixed space ->", repr(e))
report("ASMLinesmoothPC codims 0,1 on V (u only)", patches(asm.ASMLinesmoothPC, V, "lu_", codims="0,1"), f"cells {6*nz+2}, vertices {2*nz}")
report("ASMLinesmoothPC codims 1 on V (u only)", patches(asm.ASMLinesmoothPC, V, "lu1_", codims="1"), f"{2*nz}")
# 'reduced star' = dofs on the base vertex only
Vidx = asm.get_local_ises_indices(Wm)
dm = mesh.topology_dm
vstart, vend = dm.getDepthStratum(0)
red = [np.array(asm.get_entity_dofs(Wm, Vidx, [p])) for p in range(vstart, vend)]
report("get_entity_dofs(vertex) (mixed) [reduced star]", red, f"{2*nz}")

# ---- divergence-free content of each patch (u-dofs only, so enumerate on V) ----
Q = FunctionSpace(mesh, 'DG', deg-1)
B = assemble(TestFunction(Q)*div(TrialFunction(V))*dx).M.handle[:, :]
print("div matrix", B.shape, " kernel dim", B.shape[1]-np.linalg.matrix_rank(B), " expected kNx(kNz-1)+1 =", deg*nx*(deg*nz-1)+1)
def kerdim(idx):
    sub = B[:, idx]
    return len(idx) - np.linalg.matrix_rank(sub)
def kreport(name, ps, expect):
    ks = sorted(set(kerdim(p) for p in ps))
    print(f"  ker(div) inside patch, {name}: {ks}   expected {expect}")
kreport("star dim0 (u)", patches(asm.ASMStarPC, V, "ks0_", construct_dim=0), f"(2k-1)(kNz-1) = {(2*deg-1)*(deg*nz-1)}")
kreport("star dim1 = cell column (u)", patches(asm.ASMStarPC, V, "ks1_", construct_dim=1), f"(k-1)(kNz-1) = {(deg-1)*(deg*nz-1)}")
kreport("vanka dim1 = single column (u)", patches(asm.ASMVankaPC, V, "kv1_", construct_dim=1), f"(k-1)(kNz-1) = {(deg-1)*(deg*nz-1)}")
kreport("facet line (u)", patches(asm.ASMLinesmoothPC, V, "kl1_", codims="1"), "0")
kreport("extruded vertex star (u)", patches(asm.ASMExtrudedStarPC, V, "ke0_", construct_dim=0), "9 interior, 3 at top/bottom")
# ---- u-dofs of the facet line do not couple across layers in the reduced form ----
m_, c_ = Constant(37.0), Constant(1e3)
u, w = TrialFunction(V), TestFunction(V)
A = assemble((u[0]*w[0] + m_*u[1]*w[1] + c_*div(u)*div(w))*dx).M.handle[:, :]
for p in patches(asm.ASMLinesmoothPC, V, "al1_", codims="1")[:1]:
    sub = A[np.ix_(p, p)]
    nzs = (np.abs(sub) > 1e-12*np.abs(sub).max()).sum()
    print("facet-line block: size", len(p), " nonzeros", nzs, " (2x2 block diagonal would give", 2*len(p), ")")
