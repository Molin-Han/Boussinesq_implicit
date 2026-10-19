"""
02_firedrake_crosscheck.py

Checks that the tensor-product numpy model (fe1d.py) used in theory_kernel.md is
the same object as what the solver script builds in Firedrake:

  * mesh built exactly as in Nonlinear_Boussinesq_Irk_SC.py
    (PeriodicIntervalMesh -> MeshHierarchy -> ExtrudedMeshHierarchy ->
     utils.high_dim_mesh_hierarchy), velocity space utils.extrude_RT(mesh, k),
    DirichletBC on "top" and "bottom";
  * bilinear form a(u,w) = (u_x,w_x) + m (u_z,w_z) + c (div u, div w);
  * patches taken from firedrake.ASMStarPC (construct_dim 0) and
    firedrake.ASMVankaPC (construct_dim 1) through their own get_patches().

Compared quantities (Firedrake vs numpy model):
  dim V_h, dim N_h, patch sizes, dim N_h∩V_j, dim sum_j N_h∩V_j,
  dim N_H + sum_j N_h∩V_j (coarse space through firedrake.TransferManager().prolong),
  and the full spectrum of the one-level and two-level additive Schwarz
  preconditioned operators (these are basis independent, so they must agree to
  round-off if and only if operator, patch subspaces and coarse space agree).

Run (serial):  python 02_firedrake_crosscheck.py > 02_firedrake_crosscheck.txt
"""
import sys
import numpy as np
sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")
from firedrake import *                                      # noqa
from firedrake.preconditioners.asm import ASMStarPC, ASMVankaPC
from firedrake.petsc import PETSc
import utils
from fe1d import Model2D, rank, nullspace, asm_inverse, gen_eigs

Print = PETSc.Sys.Print


def build(nx_base, nz, Lx, Lz, k, refinement=1):
    dist = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
    m = PeriodicIntervalMesh(nx_base, Lx, distribution_parameters=dist)
    mh = MeshHierarchy(m, refinement_levels=refinement)
    hierarchy = ExtrudedMeshHierarchy(mh, Lz, layers=[nz] * (refinement + 1), extrusion_type='uniform')
    new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
    return new_mh


def fd_operator(mesh, k, mpar, cpar):
    V = utils.extrude_RT(mesh, k=k)
    u, w = TrialFunction(V), TestFunction(V)
    kz = as_vector([0., 0., 1.])
    a = inner(u, w) * dx + Constant(mpar - 1.0) * inner(u, kz) * inner(w, kz) * dx \
        + Constant(cpar) * div(u) * div(w) * dx
    bcs = [DirichletBC(V, as_vector([0., 0., 0.]), "top"),
           DirichletBC(V, as_vector([0., 0., 0.]), "bottom")]
    A = assemble(a, bcs=bcs, mat_type="aij").petscmat
    A = A.convert("dense").getDenseArray().copy()
    Q = FunctionSpace(mesh, "DG", k - 1)
    q = TestFunction(Q)
    B = assemble(q * div(u) * dx, mat_type="aij").petscmat.convert("dense").getDenseArray().copy()
    bcn = np.unique(np.concatenate([bc.nodes for bc in bcs]))
    free = np.setdiff1d(np.arange(V.dim()), bcn)
    return V, A, B, free, bcn


def fd_patches(V, kind):
    if kind == "star":
        pc = ASMStarPC.__new__(ASMStarPC)
        pc.prefix = "crosscheck_star_"
        PETSc.Options()["crosscheck_star_construct_dim"] = 0
    else:
        pc = ASMVankaPC.__new__(ASMVankaPC)
        pc.prefix = "crosscheck_vanka_"
        PETSc.Options()["crosscheck_vanka_construct_dim"] = 1
    ises = pc.get_patches(V)
    return [np.array(i.indices, dtype=int) for i in ises]


def fd_prolongation(Vc, Vf):
    """matrix of the prolongation Vc -> Vf that PCMG uses in the solver, column by
    column.  The solver passes no transfer manager, so Firedrake uses the default
    firedrake.TransferManager(); for this element (EnrichedElement, Piola mapped,
    immersed mesh) that is the 'embedded' transfer: project to broken DG, prolong
    natively, project back cellwise and average."""
    tm = TransferManager()
    P = np.zeros((Vf.dim(), Vc.dim()))
    uc, uf = Function(Vc), Function(Vf)
    for j in range(Vc.dim()):
        uc.dat.data[:] = 0.0
        uc.dat.data[j] = 1.0
        tm.prolong(uc, uf)
        P[:, j] = uf.dat.data_ro
    return P


def kernel_in(B, idx, n):
    Nj = nullspace(B[:, idx])
    out = np.zeros((n, Nj.shape[1]))
    out[idx, :] = Nj
    return out


def main():
    Lx, Lz = 300e3, 4e3
    dt, N = 50.0, 1e-2
    mpar = 1 + (dt * N) ** 2 / 4
    CM2 = max(Lx ** 2 / (4 * np.pi ** 2), mpar * Lz ** 2 / np.pi ** 2)
    allok = True
    for k in (1, 2):
        for Kfac in (1.0, 100.0):
            cpar = CM2 / Kfac
            nx_base, nz = 4, 3
            mh = build(nx_base, nz, Lx, Lz, k, refinement=1)
            Vf, A, B, free, bcn = fd_operator(mh[-1], k, mpar, cpar)
            Vc, Ac, Bc, freec, bcnc = fd_operator(mh[0], k, mpar, cpar)
            Nx = 2 * nx_base
            nfree = len(free)
            # restrict to free dofs
            pos = -np.ones(Vf.dim(), dtype=int)
            pos[free] = np.arange(nfree)
            Af = A[np.ix_(free, free)]
            Bf = B[:, free]
            P = fd_prolongation(Vc, Vf)[np.ix_(free, freec)]
            # numpy model
            fine = Model2D(Nx, nz, Lx, Lz, k)
            coarse = Model2D(Nx // 2, nz, Lx, Lz, k)
            Am = fine.operator(mpar, cpar)
            Pm = fine.prolongation_from(coarse)
            Bm = fine.div_matrix()
            Print(f"\n=== k = {k}, K = {Kfac:g}, Nx = {Nx}, nz = {nz}, c/dx^2 = {cpar/(Lx/Nx)**2:.4g} ===")
            Print(f"dim V_h (free dofs): firedrake {nfree}, model {fine.n}; Dirichlet dofs removed: {len(bcn)}")
            ev_f = np.linalg.eigvalsh(Af)
            ev_m = np.linalg.eigvalsh(Am)
            # the two matrices are similar only up to the basis; compare the generalised
            # spectrum A x = lambda M x with the mass matrix instead (basis independent)
            Vfun = utils.extrude_RT(mh[-1], k=k)
            uu, ww = TrialFunction(Vfun), TestFunction(Vfun)
            Mfd = assemble(inner(uu, ww) * dx, mat_type="aij").petscmat.convert("dense").getDenseArray()
            Mfd = Mfd[np.ix_(free, free)]
            import scipy.linalg as sla
            g_f = sla.eigh(Af, Mfd, eigvals_only=True)
            g_m = sla.eigh(Am, fine.mass(1.0), eigvals_only=True)
            err = np.max(np.abs(g_f - g_m) / np.abs(g_m))
            Print(f"generalised spectrum of (A, mass): max rel. difference firedrake/model = {err:.2e}")
            allok &= err < 1e-8
            dimN_f = nullspace(Bf).shape[1]
            dimN_m = nullspace(Bm).shape[1]
            Print(f"dim N_h: firedrake {dimN_f}, model {dimN_m}, formula {k*Nx*(k*nz-1)+1}")
            allok &= dimN_f == dimN_m
            NH_f = P @ nullspace(Bc[:, freec])
            Print(f"coarse kernel prolonged by TransferManager.prolong: max |div| = {np.abs(Bf @ NH_f).max():.2e} "
                  f"(relative to max|B| = {np.abs(Bf).max():.2e})")
            for kind in ("star", "vanka"):
                pats_all = fd_patches(Vf, kind)
                sizes_all = sorted(set(len(p) for p in pats_all))
                pats = [pos[p][pos[p] >= 0] for p in pats_all]
                sizes = sorted(set(len(p) for p in pats))
                pm = fine.patches(kind)
                Print(f"[{kind}] number of patches: firedrake {len(pats)}, model {len(pm)}; "
                      f"patch size incl. Dirichlet dofs {sizes_all}, free dofs {sizes}, model {sorted(set(len(p) for p in pm))}")
                Kj = [kernel_in(Bf, idx, nfree) for idx in pats]
                dNj = sorted(set(K_.shape[1] for K_ in Kj))
                allK = np.hstack(Kj)
                s1 = rank(allK) if allK.shape[1] else 0
                s2 = rank(np.hstack([NH_f, allK]))
                Print(f"[{kind}] dim N_h∩V_j = {dNj}, dim sum_j N_h∩V_j = {s1}, "
                      f"dim N_H + sum_j N_h∩V_j = {s2}, dim N_h = {dimN_f}")
                # one-level and two-level additive Schwarz spectra
                e1f = gen_eigs(asm_inverse(Af, pats), Af)
                e1m = gen_eigs(asm_inverse(Am, pm), Am)
                AH_f = P.T @ Af @ P
                AH_m = Pm.T @ Am @ Pm
                e2f = gen_eigs(asm_inverse(Af, pats) + P @ np.linalg.solve(AH_f, P.T), Af)
                e2m = gen_eigs(asm_inverse(Am, pm) + Pm @ np.linalg.solve(AH_m, Pm.T), Am)
                d1 = np.max(np.abs(e1f - e1m) / np.abs(e1m))
                d2 = np.max(np.abs(e2f - e2m) / np.abs(e2m))
                Print(f"[{kind}] one-level ASM: lambda_min = {e1f[0]:.6e} (model {e1m[0]:.6e}), "
                      f"lambda_max = {e1f[-1]:.6f} (model {e1m[-1]:.6f}), max rel. diff of all eigenvalues {d1:.1e}")
                Print(f"[{kind}] two-level ASM: lambda_min = {e2f[0]:.6e} (model {e2m[0]:.6e}), "
                      f"lambda_max = {e2f[-1]:.6f} (model {e2m[-1]:.6f}), max rel. diff of all eigenvalues {d2:.1e}")
                allok &= (d1 < 1e-6) and (d2 < 1e-6)
            # is the rediscretised coarse operator the Galerkin one?
            dG = np.max(np.abs(AH_f - Ac[np.ix_(freec, freec)])) / np.max(np.abs(AH_f))
            Print(f"coarse operator: |P^T A P - A_H(rediscretised)| / |P^T A P| = {dG:.2e}")
    Print(f"\nFiredrake and the numpy model agree in every comparison: {allok}")


if __name__ == "__main__":
    main()
