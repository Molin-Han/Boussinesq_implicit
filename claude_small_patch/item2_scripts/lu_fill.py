# Throwaway: LU fill (RCM ordering, as the solver's sub_pc_factor_mat_ordering_type rcm) of each candidate
# patch of the (u, u_y, b) block at Nz = 20, plus the DG(z) nodal variant behind the diagonal facet block.
import sys
sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla
from scipy.sparse.csgraph import reverse_cuthill_mckee
from firedrake import *
from firedrake.petsc import PETSc
from firedrake.preconditioners import asm
import utils

# DG_1 on an interval: is the mass matrix diagonal (Gauss-Legendre nodes)?
im = IntervalMesh(1, 1.0)
D = FunctionSpace(im, "DG", 1)
Md = assemble(TestFunction(D)*TrialFunction(D)*dx).M.handle[:, :]
print("DG1(interval) mass matrix:\n", np.round(Md, 6), "\n variant:", D.finat_element.fiat_equivalent.__class__.__name__,
      getattr(D.finat_element, "variant", None))

nx, nz, L, H, deg = 4, 20, 15000.0, 4000.0, 2       # dx = 3750 m, dz = 200 m (nominal cell)
base = PeriodicIntervalMesh(nx, L)
mesh = ExtrudedMesh(base, nz, layer_height=H/nz)
V = utils.extrude_RT(mesh, k=deg); Vy = FunctionSpace(mesh, 'DG', deg-1); Vb = utils.W_theta(mesh, k=deg)
Wm = V * Vy * Vb
dt, N2, c = Constant(50.0), Constant(1e-4), Constant(1e3)
(u, uy, b), (w, wy, q) = TrialFunctions(Wm), TestFunctions(Wm)
kk = as_vector([0., 1.])
a = (inner(u, w) + c*div(u)*div(w) - dt/2*inner(w, kk)*b + q*b + dt*N2/2*q*inner(u, kk) + uy*wy)*dx
bcs = [DirichletBC(Wm.sub(0), as_vector([0., 0.]), "top"), DirichletBC(Wm.sub(0), as_vector([0., 0.]), "bottom")]
A = assemble(a, bcs=bcs).M.handle
ai, aj, av = A.getValuesCSR()
A = sp.csr_matrix((av, aj, ai), shape=A.getSize())
print("block size", A.shape, "nnz", A.nnz)

opts = PETSc.Options()
def patches(cls, prefix, **kw):
    for k_, v_ in kw.items(): opts[prefix + k_] = v_
    obj = cls.__new__(cls); obj.prefix = prefix
    out = [np.array(i.indices) for i in obj.get_patches(Wm)]
    for k_ in kw: opts.delValue(prefix + k_)
    return out

def fill(p):
    S = A[p][:, p].tocsc()
    perm = reverse_cuthill_mckee(S, symmetric_mode=True)
    S = S[perm][:, perm].tocsc()
    lu = spla.splu(S, permc_spec="NATURAL", diag_pivot_thresh=0.0, options=dict(SymmetricMode=True))
    return lu.L.nnz + lu.U.nnz - len(p)   # count the unit diagonal of L once

rows = [("ASMStarPC dim0  (committed star, 2 columns + facet line)", patches(asm.ASMStarPC, "s0_", construct_dim=0)),
        ("ASMVankaPC dim1 (single column incl. both facet lines)", patches(asm.ASMVankaPC, "v1_", construct_dim=1)),
        ("ASMStarPC dim1  (cell column, no facet lines)", patches(asm.ASMStarPC, "s1_", construct_dim=1)),
        ("facet line (dofs on base vertex)", patches(asm.ASMVankaPC, "v0_", construct_dim=0, exclude_subspaces="0,1,2", include_type="entity")),
        ("ASMExtrudedStarPC dim0 (2D vertex stars)", patches(asm.ASMExtrudedStarPC, "e0_", construct_dim=0))]
print(f"{'patch family':60s} {'#patches':>8s} {'sizes':>12s} {'LU nnz/patch':>14s} {'LU nnz per base column':>22s}")
for name, ps in rows:
    fills = [fill(p) for p in ps]
    sizes = sorted(set(len(p) for p in ps))
    print(f"{name:60s} {len(ps):8d} {str(sizes):>12s} {str(sorted(set(fills))):>14s} {sum(fills)/nx:22.0f}")
# facet-line block structure
p = rows[3][1][0]; S = A[p][:, p]
print("facet-line block: size", len(p), " nnz", S.nnz, " max |offdiag| =", np.abs(S - sp.diags(S.diagonal())).max())
# Cheap 1-step cost proxy: nnz(L+U) per column is proportional to the triangular-solve flops of one patch application.
