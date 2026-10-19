'''
Common helpers for the algebraic verification of the column (vanka) patch against the
star patch on the ACTUAL discretisation of Nonlinear_Boussinesq_Irk_SC.py.

Everything here runs in serial on small meshes and works with dense/sparse
scipy matrices extracted from Firedrake.

What is built exactly as in the solver script
    * mesh hierarchy: PeriodicIntervalMesh -> MeshHierarchy -> ExtrudedMeshHierarchy
      -> utils.high_dim_mesh_hierarchy(dim=3)                       (build_hierarchy)
    * spaces: utils.extrude_RT, DG_{k-1}, utils.W_theta, DG_{k-1}  (build_spaces)
    * patches: firedrake.ASMStarPC / firedrake.ASMVankaPC .get_patches (patch_index_sets)
    * shift: claude_paper_shift.shift_value(K, dt, L_x, L_z)       (imported from the repo)

What is written by hand here (and checked against the solver objects by
01_verify_operators_against_solver.py)
    * the field-1 bilinear form linearised about rest with f = 0   (field1_form)
    * the weighted H(div) velocity form after eliminating b        (velocity_form)
    * the divergence form                                          (divergence_form)
'''
import os
import sys

REPO = "/home/mh2520/molin/phd/yr2/Boussinesq_implicit"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import warnings
warnings.filterwarnings("ignore", message="Discontinuous Lagrange element requested")

import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla

from firedrake import *                      # noqa: F401,F403
from firedrake.petsc import PETSc
from firedrake.preconditioners.asm import ASMStarPC, ASMVankaPC
import firedrake.logging as fd_logging

import utils                                 # the repository's utils.py
import claude_paper_shift as cps             # the repository's shift formula (single source of truth)

N_FREQ = cps.N_FREQ                          # 1e-2 1/s, = sqrt(utils.buo_freq())
NOMINAL_LENGTH = cps.LENGTH                  # 3e5 m
NOMINAL_HEIGHT = cps.HEIGHT                  # 4e3 m


# ----------------------------------------------------------------------------------------
# parameters
# ----------------------------------------------------------------------------------------
def m_factor(dt):
    '''m = 1 + (dt N)^2 / 4.'''
    return cps.m_factor(dt)


def shift_value(K, dt, length, height):
    '''delta = K * delta_thr, rounded exactly as the command generator of the sweeps does.'''
    return cps.shift_value(K, dt, length, height)


def c_factor(dt, delta):
    '''Weight of the (div u, div w) term of the velocity operator, c = dt / (2 delta).'''
    return dt / (2.0 * delta)


# ----------------------------------------------------------------------------------------
# mesh and spaces, as in the solver script
# ----------------------------------------------------------------------------------------
def build_hierarchy(nx, nz, length, height, refinement=2):
    '''
    Lines 79-83 of Nonlinear_Boussinesq_Irk_SC.py.  nx is the number of cells of the
    COARSEST base mesh; the finest level has nx * 2**refinement columns.
    '''
    distribution_parameters = {"partition": True,
                               "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
    m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
    mh = MeshHierarchy(m, refinement_levels=refinement)
    hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (refinement + 1),
                                      extrusion_type='uniform')
    new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
    return new_mh


def build_spaces(mesh, deg):
    '''Lines 101-105 of the solver script.  Returns (V_2D, Vy, Vb, Pressure, W, W1) with
    W = V_2D * Vy * Vb * Pressure (the solver space) and W1 = V_2D * Vy * Vb (field 1).'''
    V_2D = utils.extrude_RT(mesh, k=deg)
    Vy = FunctionSpace(mesh, 'DG', deg - 1)
    Pressure = FunctionSpace(mesh, 'DG', deg - 1)
    Vb = utils.W_theta(mesh, k=deg)
    W = V_2D * Vy * Vb * Pressure
    W1 = V_2D * Vy * Vb
    return V_2D, Vy, Vb, Pressure, W, W1


def velocity_bcs(V):
    '''u.n = 0 on top and bottom, as lines 125-127 / 165-167 of the solver script.'''
    zero = as_vector([0., 0., 0.])
    return [DirichletBC(V, zero, "top"), DirichletBC(V, zero, "bottom")]


# ----------------------------------------------------------------------------------------
# hand-written forms (linearised about rest, f = 0)
# ----------------------------------------------------------------------------------------
def field1_form(W1, dt, delta, N2=N_FREQ ** 2, div_of_full_u=True):
    '''
    Field-1 block (u, u_y, b) of the auxiliary operator of HDivSchurPC for the implicit
    midpoint rule (unknown = stage derivative k, state = state0 + dt/2 k), linearised
    about the state of rest with f = 0:

        (w, u) + (w_y, u_y) + (dt/2)(1/delta)(div w, div u) - (dt/2)(w.k, b)
        + (q, b) + (dt/2) N^2 (q, k.u)

    div_of_full_u=True uses div(u_xz + u_y j) as the script does (u_y j has zero surface
    divergence on the immersed slice); False uses div(u_xz) only.
    '''
    u_xz, u_y, b = TrialFunctions(W1)
    w_xz, w_y, q = TestFunctions(W1)
    kvec = as_vector([0., 0., 1.])
    jvec = as_vector([0., 1., 0.])
    u = u_xz + u_y * jvec
    w = w_xz + w_y * jvec
    half = Constant(dt / 2.0)
    shift = Constant(delta)
    if div_of_full_u:
        divu, divw = div(u), div(w)
    else:
        divu, divw = div(u_xz), div(w_xz)
    a = (inner(w, u) * dx
         + half / shift * divw * divu * dx
         - half * inner(w, kvec) * b * dx
         + q * b * dx
         + half * Constant(N2) * q * inner(kvec, u) * dx)
    return a, velocity_bcs(W1.sub(0))


def velocity_form(V, m, c):
    '''a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w) on the in-plane velocity space.'''
    u = TrialFunction(V)
    w = TestFunction(V)
    ivec = as_vector([1., 0., 0.])
    kvec = as_vector([0., 0., 1.])
    a = (inner(u, ivec) * inner(w, ivec) * dx
         + Constant(m) * inner(u, kvec) * inner(w, kvec) * dx
         + Constant(c) * div(u) * div(w) * dx)
    return a, velocity_bcs(V)


def divergence_form(V, Q):
    '''B[i, j] = (phi_i, div u_j), phi in the pressure space Q, u in V.'''
    u = TrialFunction(V)
    phi = TestFunction(Q)
    return phi * div(u) * dx


# ----------------------------------------------------------------------------------------
# matrices
# ----------------------------------------------------------------------------------------
def petsc_to_scipy(mat):
    indptr, indices, data = mat.getValuesCSR()
    return sp.csr_matrix((data, indices, indptr), shape=mat.getSize())


def assemble_scipy(a, bcs=None):
    A = assemble(a, bcs=bcs, mat_type='aij')
    return petsc_to_scipy(A.petscmat)


def field_offsets(W):
    '''Serial global index ranges of the fields of a mixed space, as arrays of indices.'''
    return [np.asarray(iset.indices).copy() for iset in W.dof_dset.field_ises]


def bc_dofs_velocity(V):
    '''Indices (in the numbering of V) of the DoFs constrained by u.n = 0 on top/bottom.'''
    nodes = np.concatenate([bc.nodes for bc in velocity_bcs(V)])
    return np.unique(nodes)


# ----------------------------------------------------------------------------------------
# patches from the real Firedrake classes
# ----------------------------------------------------------------------------------------
class _quiet_firedrake:
    '''Silence the "applying ASM...PC on an extruded mesh" warnings.'''
    def __enter__(self):
        self._lvl = fd_logging.logger.level
        fd_logging.set_log_level(fd_logging.ERROR)

    def __exit__(self, *args):
        fd_logging.set_log_level(self._lvl)


def patch_index_sets(V, kind):
    '''
    Patch index sets of the space V (mixed or not) in the numbering of the assembled
    monolithic matrix, obtained by calling the get_patches method of the real classes

        kind = "star"  : firedrake.ASMStarPC  with pc_star_construct_dim 0
        kind = "vanka" : firedrake.ASMVankaPC with pc_vanka_construct_dim 1

    followed by the local-to-global map exactly as ASMPatchPC.initialize does.
    Returns a list of integer arrays (one per patch, in the order Firedrake builds them).
    '''
    opts = PETSc.Options()
    if kind == "star":
        pc = ASMStarPC()
        pc.prefix = "fa_star_"
        opts["fa_star_construct_dim"] = 0
    elif kind == "vanka":
        pc = ASMVankaPC()
        pc.prefix = "fa_vanka_"
        opts["fa_vanka_construct_dim"] = 1
    else:
        raise ValueError(kind)
    with _quiet_firedrake():
        ises = pc.get_patches(V)
    lgmap = V.dof_dset.lgmap
    out = []
    for iset in ises:
        g = lgmap.applyIS(iset)
        out.append(np.asarray(g.indices).copy())
    return out


def restrict_patches(patches, keep):
    '''
    Restrict patches given in a global numbering to the DoFs listed in `keep` and
    renumber them 0..len(keep)-1.  Used to pass from (u, u_y, b) patches to the patches
    of the free velocity DoFs.
    '''
    keep = np.asarray(keep)
    lookup = -np.ones(int(max(keep.max(), max(p.max() for p in patches))) + 1, dtype=np.int64)
    lookup[keep] = np.arange(len(keep))
    out = []
    for p in patches:
        q = lookup[p]
        out.append(q[q >= 0])
    return out


# ----------------------------------------------------------------------------------------
# additive Schwarz and spectra
# ----------------------------------------------------------------------------------------
def asm_inverse_dense(A, patches):
    '''Dense matrix of sum_j R_j^T (R_j A R_j^T)^{-1} R_j for a dense or sparse A.'''
    n = A.shape[0]
    Ad = A.toarray() if sp.issparse(A) else np.asarray(A)
    Pinv = np.zeros((n, n))
    for p in patches:
        if len(p) == 0:
            continue
        Aj = Ad[np.ix_(p, p)]
        Pinv[np.ix_(p, p)] += np.linalg.inv(Aj)
    return Pinv


def spd_preconditioned_eigs(S, Pinv, return_vectors=False):
    '''
    Eigenvalues of Pinv S for SPD S and SPD Pinv, computed from the symmetric matrix
    L^T Pinv L with S = L L^T.  If return_vectors, also returns the eigenvectors x of
    Pinv S x = lambda x, normalised so that x^T S x = 1.
    '''
    Sd = S.toarray() if sp.issparse(S) else np.asarray(S)
    Sd = 0.5 * (Sd + Sd.T)
    L = np.linalg.cholesky(Sd)
    T = L.T @ Pinv @ L
    T = 0.5 * (T + T.T)
    if not return_vectors:
        return np.linalg.eigvalsh(T)
    lam, Y = np.linalg.eigh(T)
    # Pinv S x = lam x with x = L^{-T} y ; then x^T S x = y^T y = 1
    X = sla.solve_triangular(L.T, Y, lower=False)
    return lam, X


def spectrum_summary(lam, thresholds=(0.5, 0.1, 0.01, 0.001)):
    lam = np.sort(np.asarray(lam))
    out = {"n": len(lam), "lmin": lam[0], "lmax": lam[-1], "cond": lam[-1] / lam[0]}
    for th in thresholds:
        out[f"n_lt_{th:g}"] = int(np.sum(lam < th))
    return out


def null_space_dim(B, rtol=1e-10):
    '''Dimension of the kernel of a dense matrix B (n_rows x n_cols) and its singular values.'''
    if B.shape[1] == 0:
        return 0, np.zeros(0)
    s = np.linalg.svd(B, compute_uv=False)
    smax = s[0] if len(s) and s[0] > 0 else 1.0
    rank = int(np.sum(s > rtol * smax))
    return B.shape[1] - rank, s


def null_space_basis(B, rtol=1e-10):
    '''Orthonormal basis (columns) of the kernel of the dense matrix B.'''
    if B.shape[1] == 0:
        return np.zeros((0, 0))
    U, s, Vt = np.linalg.svd(B, full_matrices=True)
    smax = s[0] if len(s) and s[0] > 0 else 1.0
    rank = int(np.sum(s > rtol * smax))
    return Vt[rank:].T.copy()


def write_table(path, header, rows, fmt=None):
    '''Write a simple CSV file.'''
    with open(path, "w") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(_fmt(x) for x in r) + "\n")


def _fmt(x):
    if isinstance(x, (float, np.floating)):
        return f"{x:.6e}"
    return str(x)


# ----------------------------------------------------------------------------------------
# the weighted H(div) velocity problem on one mesh, assembled once in pieces
# ----------------------------------------------------------------------------------------
class VelocityProblem:
    '''
    Pieces of the velocity operator of field 1 (after eliminating b, rest state, f = 0)

        S(m, c) = M_x + m M_z + c D,
        M_x = (u_x, w_x),  M_z = (u_z, w_z),  D = (div u, div w),

    on the finest mesh of the hierarchy built as in the solver script, restricted to the
    free velocity DoFs (u.n = 0 on top and bottom removed), together with the velocity part
    of the star and vanka patches of the mixed space (u, u_y, b) built by the real
    Firedrake classes.
    '''
    def __init__(self, nx, nz, length, height, deg, refinement=2, level=-1):
        self.nx, self.nz, self.length, self.height, self.deg = nx, nz, length, height, deg
        self.refinement = refinement
        self.mh = build_hierarchy(nx, nz, length, height, refinement=refinement)
        self.mesh = self.mh[level]
        V, Vy, Vb, Q, W, W1 = build_spaces(self.mesh, deg)
        self.V, self.Vy, self.Vb, self.Q, self.W, self.W1 = V, Vy, Vb, Q, W, W1
        self.ncol = self.mesh.num_cells()
        self.dx = length / self.ncol
        self.dz = height / nz
        u = TrialFunction(V)
        w = TestFunction(V)
        ivec = as_vector([1., 0., 0.])
        kvec = as_vector([0., 0., 1.])
        bcn = bc_dofs_velocity(V)
        self.bc_nodes = bcn
        self.free = np.setdiff1d(np.arange(V.dim()), bcn)
        fr = self.free
        self.Mx = assemble_scipy(inner(u, ivec) * inner(w, ivec) * dx)[fr][:, fr].tocsr()
        self.Mz = assemble_scipy(inner(u, kvec) * inner(w, kvec) * dx)[fr][:, fr].tocsr()
        self.D = assemble_scipy(div(u) * div(w) * dx)[fr][:, fr].tocsr()
        self.Dzz = assemble_scipy(grad(u)[2, 2] * grad(w)[2, 2] * dx)[fr][:, fr].tocsr()
        self.B = assemble_scipy(divergence_form(V, Q))[:, fr].tocsr()
        self.Mp = assemble_scipy(TrialFunction(Q) * TestFunction(Q) * dx)
        self.n = len(fr)
        f1 = field_offsets(W1)
        assert np.array_equal(f1[0], np.arange(V.dim()))
        self.patches_mixed = {kind: patch_index_sets(W1, kind) for kind in ("star", "vanka")}
        self.patches = {kind: restrict_patches(self.patches_mixed[kind], fr)
                        for kind in ("star", "vanka")}

    def S(self, m, c):
        return (self.Mx + m * self.Mz + c * self.D).toarray()

    def S_degraded(self, m, c):
        '''The operator the hypothesis has in mind: div u replaced by d_z u_z.'''
        return (self.Mx + m * self.Mz + c * self.Dzz).toarray()

    def params(self, dt, K, length_for_shift=None, height_for_shift=None):
        '''(delta, m, c) of the experiments: delta = K delta_thr(dt, L_x, L_z).'''
        L = self.length if length_for_shift is None else length_for_shift
        H = self.height if height_for_shift is None else height_for_shift
        delta = shift_value(K, dt, L, H)
        return delta, m_factor(dt), c_factor(dt, delta)
