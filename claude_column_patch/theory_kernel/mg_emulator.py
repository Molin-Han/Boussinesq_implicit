"""
mg_emulator.py -- numpy/scipy emulation of the field-1 multigrid of the solver,
applied to the symmetric velocity operator

    a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w)

(buoyancy eliminated, f = 0, linearised about rest), discretised exactly as in
the solver (tensor-product model of fe1d.py, checked against Firedrake in
02_firedrake_crosscheck.py).

Algorithm (the PETSc options of the solver script):
    pc_mg_type full, pc_mg_cycle_type v, refinement+1 levels (coarsening in x only),
    mg_levels: GMRES (left preconditioned, 6 iterations, nonzero initial guess),
               preconditioned by the additive Schwarz method with star or
               single-column patches, exact patch solves,
    mg_coarse: exact solve,
    restriction = transpose of the natural injection, coarse operators = rediscretised
    (= Galerkin, the spaces are nested).
Outer iteration: FGMRES(30), right preconditioned by ONE multigrid application.

This is an emulation of the multigrid component only.  In the solver the same
multigrid acts on the nonsymmetric block (u, u_y, b) inside a Schur-complement
factorisation; iteration counts of the two are not expected to be equal.
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import scipy.sparse.linalg as spla
from fe1d import Model2D, prolongation_x


def sparse_operator(mod, m, c):
    X, Z = mod.X, mod.Z
    S = sp.csr_matrix
    M = sp.block_diag([sp.kron(S(X.Mc), S(Z.Md)), m * sp.kron(S(X.Md), S(Z.Mc))], format="csr")
    B = sp.hstack([sp.kron(S(X.D), S(Z.Md)), sp.kron(S(X.Md), S(Z.D))], format="csr")
    MQ = sp.kron(S(X.Md), S(Z.Md)).diagonal()
    A = M + c * (B.T @ sp.diags(1.0 / MQ) @ B)
    A = A.tocsr()
    A.eliminate_zeros()
    return A, M, B, MQ


def sparse_prolongation(coarse, fine):
    Pc, Pd = prolongation_x(coarse.X, fine.X)
    S = sp.csr_matrix
    return sp.block_diag([sp.kron(S(Pc), sp.identity(fine.Z.nd)),
                          sp.kron(S(Pd), sp.identity(fine.Z.nc))], format="csr")


class ASM:
    def __init__(self, A, patches):
        self.patches = patches
        self.lu = []
        Ac = A.tocsc()
        for idx in patches:
            Aj = A[idx][:, idx].toarray()
            self.lu.append(sla.lu_factor(Aj))
        self.n = A.shape[0]

    def apply(self, r):
        y = np.zeros(self.n)
        for idx, lu in zip(self.patches, self.lu):
            y[idx] += sla.lu_solve(lu, r[idx])
        return y


def gmres_fixed(A, Minv, b, x0, nit):
    """nit iterations of left-preconditioned GMRES for A x = b starting from x0
    (minimises |Minv (b - A x)| over x0 + Krylov space), no convergence test."""
    r = Minv(b - A @ x0)
    beta = np.linalg.norm(r)
    if beta == 0.0:
        return x0
    n = len(b)
    V = np.zeros((n, nit + 1))
    H = np.zeros((nit + 1, nit))
    V[:, 0] = r / beta
    kdim = nit
    for j in range(nit):
        w = Minv(A @ V[:, j])
        for i in range(j + 1):
            H[i, j] = V[:, i] @ w
            w -= H[i, j] * V[:, i]
        for i in range(j + 1):                       # reorthogonalise once
            t = V[:, i] @ w
            H[i, j] += t
            w -= t * V[:, i]
        H[j + 1, j] = np.linalg.norm(w)
        if H[j + 1, j] < 1e-14 * beta:
            kdim = j + 1
            break
        V[:, j + 1] = w / H[j + 1, j]
    e1 = np.zeros(kdim + 1)
    e1[0] = beta
    y = np.linalg.lstsq(H[:kdim + 1, :kdim], e1, rcond=None)[0]
    return x0 + V[:, :kdim] @ y


class Multigrid:
    def __init__(self, k, nx_coarse, nz, Lx, Lz, m, c, kind, refinement=2, smooth_its=6,
                 smoother="gmres", omega=None):
        self.levels = []
        self.kind = kind
        self.smooth_its = smooth_its
        self.smoother = smoother
        for l in range(refinement + 1):
            mod = Model2D(nx_coarse * 2 ** l, nz, Lx, Lz, k)
            A, M, B, MQ = sparse_operator(mod, m, c)
            lev = dict(mod=mod, A=A, M=M)
            if l == 0:
                lev["solve"] = spla.splu(A.tocsc()).solve
            else:
                lev["P"] = sparse_prolongation(self.levels[-1]["mod"], mod)
                lev["asm"] = ASM(A, mod.patches(kind))
            self.levels.append(lev)
        self.L = refinement
        self.A = self.levels[-1]["A"]
        self.n = self.A.shape[0]
        self.omega = omega
        self.nasm = 0

    def smooth(self, l, b, x):
        lev = self.levels[l]
        if self.smoother == "gmres":
            return gmres_fixed(lev["A"], lev["asm"].apply, b, x, self.smooth_its)
        elif self.smoother == "richardson":
            for _ in range(self.smooth_its):
                x = x + self.omega * lev["asm"].apply(b - lev["A"] @ x)
            return x
        raise ValueError

    def vcycle(self, l, b, x):
        if l == 0:
            return self.levels[0]["solve"](b)
        lev = self.levels[l]
        x = self.smooth(l, b, x)
        r = b - lev["A"] @ x
        xc = self.vcycle(l - 1, lev["P"].T @ r, np.zeros(lev["P"].shape[1]))
        x = x + lev["P"] @ xc
        x = self.smooth(l, b, x)
        return x

    def fcycle(self, b):
        bs = [None] * (self.L + 1)
        bs[self.L] = b
        for l in range(self.L, 0, -1):
            bs[l - 1] = self.levels[l]["P"].T @ bs[l]
        x = self.levels[0]["solve"](bs[0])
        for l in range(1, self.L + 1):
            x = self.levels[l]["P"] @ x
            x = self.vcycle(l, bs[l], x)
        return x

    def vcycle_top(self, b):
        return self.vcycle(self.L, b, np.zeros(self.n))


def fgmres(A, prec, b, rtol=1e-6, maxit=150, restart=30):
    """right-preconditioned flexible GMRES, zero initial guess.
    Returns (x, list of relative residual norms)."""
    n = len(b)
    x = np.zeros(n)
    bnorm = np.linalg.norm(b)
    hist = [1.0]
    its = 0
    r = b.copy()
    while its < maxit:
        beta = np.linalg.norm(r)
        V = np.zeros((n, restart + 1))
        Zs = np.zeros((n, restart))
        H = np.zeros((restart + 1, restart))
        V[:, 0] = r / beta
        kk = 0
        done = False
        for j in range(restart):
            Zs[:, j] = prec(V[:, j])
            w = A @ Zs[:, j]
            for i in range(j + 1):
                H[i, j] = V[:, i] @ w
                w -= H[i, j] * V[:, i]
            for i in range(j + 1):
                t = V[:, i] @ w
                H[i, j] += t
                w -= t * V[:, i]
            H[j + 1, j] = np.linalg.norm(w)
            V[:, j + 1] = w / H[j + 1, j] if H[j + 1, j] > 0 else 0.0
            kk = j + 1
            its += 1
            e1 = np.zeros(kk + 1)
            e1[0] = beta
            y, res, rk, sv = np.linalg.lstsq(H[:kk + 1, :kk], e1, rcond=None)
            rn = np.linalg.norm(e1 - H[:kk + 1, :kk] @ y) / bnorm
            hist.append(rn)
            if rn < rtol or its >= maxit:
                done = True
                break
        x = x + Zs[:, :kk] @ y
        r = b - A @ x
        if done:
            break
    return x, hist


def physical_rhs(mod, amp=3e-4, a=5000.0):
    """load vector of w -> (w_z, b0) with the initial buoyancy of the solver script,
    b0 = amp sin(pi z / Lz) / (1 + (x - Lx/2)^2 / a^2)   (u_x part of the load is zero)."""
    from fe1d import gauss01, lagrange_basis, legendre_basis
    X, Z, k = mod.X, mod.Z, mod.k
    xi, w = gauss01(8)
    leg = legendre_basis(k - 1)(xi)
    lag = lagrange_basis(k)[0](xi)
    fx = np.zeros(X.nd)
    for j in range(X.n):
        xq = (j + xi) * X.h
        g = 1.0 / (1.0 + (xq - mod.Lx / 2) ** 2 / a ** 2)
        fx[X.dg_cell_dofs(j)] += leg.T @ (w * X.h * g)
    fz_full = np.zeros(Z.nc_full)
    for j in range(Z.n):
        zq = (j + xi) * Z.h
        g = amp * np.sin(np.pi * zq / mod.Lz)
        fz_full[Z.cg_cell_dofs(j)] += lag.T @ (w * Z.h * g)
    fz = fz_full[Z.free]
    return np.concatenate([np.zeros(mod.n_ux), np.kron(fx, fz)])
