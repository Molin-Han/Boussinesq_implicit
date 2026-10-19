"""
fe1d.py -- building blocks for the pencil-and-paper analysis in theory_kernel.md

Tensor-product model of the velocity block of field 1 (linearised about rest,
f = 0, buoyancy eliminated):

    a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w)

on a mesh of Nx x nz rectangles dx x dz, periodic in x, u.n = 0 on top and
bottom, with

    u_x in CG_k(x) (x) DG_{k-1}(z),     u_z in DG_{k-1}(x) (x) CG_k^0(z)

which is exactly the space utils.extrude_RT(mesh, k) with the DirichletBCs of
the solver script (checked against Firedrake in 02_firedrake_crosscheck.py).

Everything is numpy/scipy only.  Conventions
  * x cells j = 0..Nx-1 lie between the vertices j and j+1 (mod Nx)
  * CG_k(x) dofs: vertex i -> i ; interior dofs of cell j -> Nx + j (k-1) + (a-1)
  * DG_{k-1}(x) dofs: cell j, local alpha -> j k + alpha   (Legendre basis)
  * 2D numbering: u_x dof (a_x, beta_z) -> a_x * n_dz + beta_z
                  u_z dof (alpha_x, b_z) -> n_ux + alpha_x * n_cz + b_z
"""
import numpy as np
from numpy.polynomial import legendre as npleg
import scipy.linalg as sla


# ----------------------------------------------------------------------------
# reference bases on [0, 1]
# ----------------------------------------------------------------------------
def lagrange_basis(k):
    """Lagrange basis of P_k on [0,1] at the equispaced nodes a/k (for k = 1, 2
    these are also the GLL nodes used by Firedrake).  Returns (val, der), each a
    function xi -> array (len(xi), k+1)."""
    nodes = np.arange(k + 1) / k
    V = np.vander(nodes, k + 1, increasing=True)      # V[a, p] = nodes_a^p
    C = np.linalg.inv(V)                               # basis_b(xi) = sum_p C[p, b] xi^p

    def val(xi):
        xi = np.atleast_1d(xi)
        P = np.vander(xi, k + 1, increasing=True)
        return P @ C

    def der(xi):
        xi = np.atleast_1d(xi)
        P = np.zeros((len(xi), k + 1))
        for p in range(1, k + 1):
            P[:, p] = p * xi ** (p - 1)
        return P @ C

    return val, der


def legendre_basis(q):
    """Legendre basis of P_q on [0,1] (orthogonal).  Returns val: xi -> (len(xi), q+1)."""
    def val(xi):
        xi = np.atleast_1d(xi)
        out = np.zeros((len(xi), q + 1))
        for n in range(q + 1):
            c = np.zeros(n + 1)
            c[n] = 1.0
            out[:, n] = npleg.legval(2 * xi - 1, c)
        return out
    return val


def gauss01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1), 0.5 * w


# ----------------------------------------------------------------------------
# 1D spaces
# ----------------------------------------------------------------------------
class Space1D:
    """CG_k and DG_{k-1} on a uniform 1D mesh with n cells of size h.

    periodic=True : CG_k periodic (dim k n)
    periodic=False: CG_k with homogeneous Dirichlet conditions at both ends
                    (dim k n - 1); the constrained vertices are removed.
    """

    def __init__(self, n, h, k, periodic):
        self.n, self.h, self.k, self.periodic = n, h, k, periodic
        self.nd = k * n                                  # dim DG_{k-1}
        if periodic:
            self.nvert = n
        else:
            self.nvert = n + 1
        self.nc_full = self.nvert + n * (k - 1)          # dim CG_k before BCs
        if periodic:
            self.free = np.arange(self.nc_full)
        else:
            self.free = np.array([i for i in range(self.nc_full) if i not in (0, n)])
        self.nc = len(self.free)
        self._assemble()

    # local -> global maps (global = numbering before removal of BC dofs)
    def cg_cell_dofs(self, j):
        k, n = self.k, self.n
        left = j
        right = (j + 1) % n if self.periodic else j + 1
        interior = [self.nvert + j * (k - 1) + (a - 1) for a in range(1, k)]
        return [left] + interior + [right]              # ordered by node a = 0..k

    def dg_cell_dofs(self, j):
        return [j * self.k + al for al in range(self.k)]

    def _assemble(self):
        k, n, h = self.k, self.n, self.h
        xi, w = gauss01(k + 2)
        val, der = lagrange_basis(k)
        leg = legendre_basis(k - 1)
        Phi, dPhi, Psi = val(xi), der(xi) / h, leg(xi)
        W = w * h
        Mc_loc = Phi.T @ (W[:, None] * Phi)
        Kc_loc = dPhi.T @ (W[:, None] * dPhi)
        Md_loc = Psi.T @ (W[:, None] * Psi)
        D_loc = Psi.T @ (W[:, None] * dPhi)             # (psi_alpha, phi_a')
        Mc = np.zeros((self.nc_full, self.nc_full))
        Kc = np.zeros_like(Mc)
        Md = np.zeros((self.nd, self.nd))
        D = np.zeros((self.nd, self.nc_full))
        for j in range(n):
            c = self.cg_cell_dofs(j)
            d = self.dg_cell_dofs(j)
            Mc[np.ix_(c, c)] += Mc_loc
            Kc[np.ix_(c, c)] += Kc_loc
            Md[np.ix_(d, d)] += Md_loc
            D[np.ix_(d, c)] += D_loc
        f = self.free
        self.Mc = Mc[np.ix_(f, f)]
        self.Kc = Kc[np.ix_(f, f)]
        self.Md = Md
        self.D = D[:, f]
        # position of a full-numbering dof in the reduced numbering (-1 if removed)
        self.full2free = -np.ones(self.nc_full, dtype=int)
        self.full2free[f] = np.arange(len(f))

    # ---- patches (x direction only; periodic) --------------------------------
    def star_patch(self, i):
        """CG and DG dofs of the star of vertex i: vertex i + interiors of the two
        adjacent cells (CG), both adjacent cells (DG)."""
        assert self.periodic
        n, k = self.n, self.k
        cells = [(i - 1) % n, i]
        cg = [i] + [self.nvert + j * (k - 1) + (a - 1) for j in cells for a in range(1, k)]
        dg = [d for j in cells for d in self.dg_cell_dofs(j)]
        return np.array(sorted(set(cg))), np.array(sorted(set(dg)))

    def vanka_patch(self, j):
        """closure of cell j: both vertices + interior (CG), the cell (DG)."""
        assert self.periodic
        cg = self.cg_cell_dofs(j)
        dg = self.dg_cell_dofs(j)
        return np.array(sorted(set(cg))), np.array(sorted(set(dg)))

    def cellint_patch(self, j):
        """interior of cell j only: interior CG dofs (none for k = 1), the cell (DG)."""
        k = self.k
        cg = [self.nvert + j * (k - 1) + (a - 1) for a in range(1, k)]
        dg = self.dg_cell_dofs(j)
        return np.array(sorted(set(cg)), dtype=int), np.array(sorted(set(dg)))


def prolongation_x(coarse, fine):
    """Natural injection CG_k(2h) -> CG_k(h) and DG_{k-1}(2h) -> DG_{k-1}(h),
    periodic, fine.n = 2 coarse.n.  Returns (Pc, Pd)."""
    assert fine.n == 2 * coarse.n and fine.k == coarse.k and fine.periodic and coarse.periodic
    k = fine.k
    val, _ = lagrange_basis(k)
    leg = legendre_basis(k - 1)
    nodes = np.arange(k + 1) / k
    Pc = np.zeros((fine.nc_full, coarse.nc_full))
    Pd = np.zeros((fine.nd, coarse.nd))
    xi, w = gauss01(k + 2)
    Psi = leg(xi)
    Md_loc = Psi.T @ (w[:, None] * Psi)
    for J in range(coarse.n):
        cJ = coarse.cg_cell_dofs(J)
        dJ = coarse.dg_cell_dofs(J)
        for s in range(2):
            j = 2 * J + s
            cj = fine.cg_cell_dofs(j)
            dj = fine.dg_cell_dofs(j)
            Xi_nodes = (s + nodes) / 2.0                 # coarse reference coordinate of the fine nodes
            Pc[np.ix_(cj, cJ)] = val(Xi_nodes)
            Xi_q = (s + xi) / 2.0
            rhs = Psi.T @ (w[:, None] * leg(Xi_q))       # (psi_fine, Psi_coarse)
            Pd[np.ix_(dj, dJ)] = np.linalg.solve(Md_loc, rhs)
    return Pc, Pd


# ----------------------------------------------------------------------------
# vertical modes: generalised SVD of d/dz : CG_k^0(z) -> DG_{k-1}(z)
# ----------------------------------------------------------------------------
def vertical_modes(nz, dz, k):
    """sigma_n (n = 1 .. k nz - 1): the singular values of d/dz from CG_k^0 to
    DG_{k-1} with respect to the L2 inner products, i.e. sigma_n^2 are the
    eigenvalues of the discrete Dirichlet Laplacian.  sigma = 0 is the extra
    (vertically constant, 'barotropic') mode of DG_{k-1}(z)."""
    Z = Space1D(nz, dz, k, periodic=False)
    lam = sla.eigh(Z.Kc, Z.Mc, eigvals_only=True)
    return np.sqrt(np.maximum(lam, 0.0))


# ----------------------------------------------------------------------------
# 1D mode operator  a_sigma((U,W),(U,W)) = |U|^2 + m |W|^2 + c |U' + sigma W|^2
# ----------------------------------------------------------------------------
def mode_operator(X, sigma, m, c):
    """Matrix of a_sigma on CG_k(x) x DG_{k-1}(x).  For sigma = 0 the W block
    does not exist (there is no u_z mode with zero vertical derivative) and the
    operator is M + c K on CG_k(x) alone."""
    if sigma == 0.0:
        return X.Mc + c * X.Kc
    nc, nd = X.nc, X.nd
    A = np.zeros((nc + nd, nc + nd))
    A[:nc, :nc] = X.Mc + c * X.Kc
    A[:nc, nc:] = c * sigma * X.D.T
    A[nc:, :nc] = c * sigma * X.D
    A[nc:, nc:] = (m + c * sigma ** 2) * X.Md
    return A


def mode_patches(X, kind, sigma):
    """index sets of the patches in the numbering of mode_operator"""
    out = []
    for p in range(X.n):
        cg, dg = getattr(X, kind + "_patch")(p)
        if sigma == 0.0:
            idx = cg
        else:
            idx = np.concatenate([cg, X.nc + dg])
        if len(idx):
            out.append(idx.astype(int))
    return out


def mode_prolongation(XH, Xh, sigma):
    Pc, Pd = prolongation_x(XH, Xh)
    if sigma == 0.0:
        return Pc
    return sla.block_diag(Pc, Pd)


def asm_inverse(A, patches):
    """dense matrix of sum_j R_j^T A_j^{-1} R_j"""
    Minv = np.zeros_like(A)
    for idx in patches:
        Minv[np.ix_(idx, idx)] += np.linalg.inv(A[np.ix_(idx, idx)])
    return Minv


def gen_eigs(Binv, A):
    """eigenvalues of Binv A for SPD A and symmetric positive semidefinite Binv,
    computed through the symmetric matrix L^T Binv L with A = L L^T."""
    L = np.linalg.cholesky(A)
    S = L.T @ Binv @ L
    S = 0.5 * (S + S.T)
    return np.linalg.eigvalsh(S)


# ----------------------------------------------------------------------------
# full 2D tensor-product operator (for the checks that do not use the mode
# decomposition)
# ----------------------------------------------------------------------------
class Model2D:
    def __init__(self, Nx, nz, Lx, Lz, k):
        self.Nx, self.nz, self.Lx, self.Lz, self.k = Nx, nz, Lx, Lz, k
        self.dx, self.dz = Lx / Nx, Lz / nz
        self.X = Space1D(Nx, self.dx, k, periodic=True)
        self.Z = Space1D(nz, self.dz, k, periodic=False)
        X, Z = self.X, self.Z
        self.n_ux = X.nc * Z.nd
        self.n_uz = X.nd * Z.nc
        self.n = self.n_ux + self.n_uz
        self.n_p = X.nd * Z.nd

    def div_matrix(self):
        """B[q, u] = (q, div u) for q in Q_h = DG(x) (x) DG(z)"""
        X, Z = self.X, self.Z
        return np.hstack([np.kron(X.D, Z.Md), np.kron(X.Md, Z.D)])

    def pressure_mass(self):
        return np.kron(self.X.Md, self.Z.Md)

    def mass(self, m=1.0):
        X, Z = self.X, self.Z
        return sla.block_diag(np.kron(X.Mc, Z.Md), m * np.kron(X.Md, Z.Mc))

    def operator(self, m, c):
        B = self.div_matrix()
        MQ = self.pressure_mass()
        return self.mass(m) + c * B.T @ np.linalg.solve(MQ, B)

    def patch(self, kind, p):
        X, Z = self.X, self.Z
        cg, dg = getattr(X, kind + "_patch")(p)
        ux = (cg[:, None] * Z.nd + np.arange(Z.nd)[None, :]).ravel()
        uz = self.n_ux + (dg[:, None] * Z.nc + np.arange(Z.nc)[None, :]).ravel()
        return np.concatenate([ux, uz]).astype(int)

    def patches(self, kind):
        return [self.patch(kind, p) for p in range(self.Nx)]

    def prolongation_from(self, coarse):
        Pc, Pd = prolongation_x(coarse.X, self.X)
        return sla.block_diag(np.kron(Pc, np.eye(self.Z.nd)), np.kron(Pd, np.eye(self.Z.nc)))

    def curl_matrix(self):
        """matrix of curl psi = (-d_z psi, d_x psi) from
        Psi_h = CG_k(x) (x) CG_k^0(z) into V_h (coefficients)."""
        X, Z = self.X, self.Z
        # d_z : CG(z) -> DG(z) and d_x : CG(x) -> DG(x) in coefficients
        Gz = np.linalg.solve(Z.Md, Z.D)
        Gx = np.linalg.solve(X.Md, X.D)
        top = -np.kron(np.eye(X.nc), Gz)
        bot = np.kron(Gx, np.eye(Z.nc))
        return np.vstack([top, bot])

    def harmonic(self):
        """coefficients of the field (1, 0)"""
        X, Z = self.X, self.Z
        one_x = np.ones(X.nc)                            # Lagrange: constant 1
        one_z = np.zeros(Z.nd)
        one_z[0::self.k] = 1.0                           # Legendre P0 coefficient
        return np.concatenate([np.kron(one_x, one_z), np.zeros(self.n_uz)])


def rank(M, tol=1e-9):
    if M.size == 0:
        return 0
    s = np.linalg.svd(M, compute_uv=False)
    if len(s) == 0 or s[0] == 0:
        return 0
    return int(np.sum(s > tol * s[0]))


def nullspace(M, tol=1e-9):
    """orthonormal basis of ker M (columns)"""
    u, s, vh = np.linalg.svd(M)
    r = int(np.sum(s > tol * s[0])) if s[0] > 0 else 0
    return vh[r:].T
