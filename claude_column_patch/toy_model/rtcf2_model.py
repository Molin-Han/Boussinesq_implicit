'''
Degree 2 version of the toy model (the element of the solver script with --degree 2):

    u_x in CG2(x) (x) DG1(z),      u_z in DG1(x) (x) CG2_0(z),

i.e. HDiv(CG2 x DG1) + HDiv(DG1 x CG2) on rectangles, periodic in x, u_z = 0 on top and
bottom, for the same bilinear form

    a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w).

1D bases:  CG2 : Lagrange nodes at the vertices and the midpoints of the cells,
           DG1 : orthonormal Legendre basis on every cell (mass matrix = identity).

Degrees of freedom attached to the entities of the base (interval) mesh:

    base vertex i : Uv[i, :]   CG2 vertex node (x)  DG1(z)          2 Nz
    base cell   i : Um[i, :]   CG2 midpoint node (x) DG1(z)         2 Nz
                    V[i, :]    DG1(x) (2 functions) (x) CG2_0(z)    2 (2 Nz - 1)

Global numbering: block i = [Uv[i], Um[i], V[i]], block size nb = 8 Nz - 2.

Patches
    star   of base vertex i : Um[i-1], V[i-1], Uv[i], Um[i], V[i]        14 Nz - 4
    column of base cell   i : Uv[i], Um[i], V[i], Uv[i+1]                10 Nz - 2
(Firedrake reports 14 Nz + 4 and 10 Nz + 2 because it also counts the constrained u_z DoFs
on the top and bottom boundary.)
'''
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla
from rt0_model import gmres_fixed, fgmres, physical_parameters

PATCHES = ('star', 'column')

MLOC = np.array([[4.0, 2.0, -1.0], [2.0, 16.0, 2.0], [-1.0, 2.0, 4.0]]) / 30.0   # times h
# derivative of the CG2 basis expressed in the orthonormal DG1 basis: (2 x 3), times 1/sqrt(h)
DLOC = np.array([[-1.0, 0.0, 1.0], [2.0 / np.sqrt(3.0), -4.0 / np.sqrt(3.0), 2.0 / np.sqrt(3.0)]])


def cg2_mass_and_derivative(n, h, periodic):
    '''CG2 mass matrix and derivative matrix D : CG2 -> DG1 (orthonormal basis) on n cells.

    periodic : nodes 2i (vertex i), 2i+1 (midpoint of cell i), 2n nodes.
    else     : nodes 0..2n, homogeneous Dirichlet at both ends, interior nodes 1..2n-1 kept
               (returned index = node - 1).
    DG1 index: 2 * cell + a.
    '''
    nn = 2 * n if periodic else 2 * n + 1
    M = sp.lil_matrix((nn, nn))
    D = sp.lil_matrix((2 * n, nn))
    for e in range(n):
        nodes = [2 * e, 2 * e + 1, (2 * e + 2) % nn if periodic else 2 * e + 2]
        for p in range(3):
            for q in range(3):
                M[nodes[p], nodes[q]] += h * MLOC[p, q]
            for a in range(2):
                D[2 * e + a, nodes[p]] += DLOC[a, p] / np.sqrt(h)
    M = M.tocsr()
    D = D.tocsr()
    if not periodic:
        keep = np.arange(1, 2 * n)
        M = M[keep][:, keep]
        D = D[:, keep]
    return M.tocsr(), D.tocsr()


def interleaved_permutation(Nx, Nz):
    '''perm[g] = field-split index of interleaved DoF g.
    field-split: U index = node_x * (2 Nz) + beta, V index = nU + (2 i + a) * (2 Nz - 1) + k.'''
    nzd = 2 * Nz
    nzc = 2 * Nz - 1
    nU = 2 * Nx * nzd
    nb = 8 * Nz - 2
    perm = np.empty(Nx * nb, dtype=np.int64)
    for i in range(Nx):
        o = i * nb
        perm[o: o + nzd] = (2 * i) * nzd + np.arange(nzd)
        perm[o + nzd: o + 2 * nzd] = (2 * i + 1) * nzd + np.arange(nzd)
        perm[o + 2 * nzd: o + nb] = nU + (2 * i) * nzc + np.arange(2 * nzc)
    return perm


def assemble(Nx, Nz, dx, dz, m, c, return_parts=False):
    Mx, Dx = cg2_mass_and_derivative(Nx, dx, periodic=True)
    Mz, Dz = cg2_mass_and_derivative(Nz, dz, periodic=False)
    Idx = sp.identity(2 * Nx)       # DG1(x) mass (orthonormal)
    Idz = sp.identity(2 * Nz)       # DG1(z) mass (orthonormal)
    M = sp.block_diag([sp.kron(Mx, Idz), sp.kron(Idx, m * Mz)]).tocsr()
    Div = sp.hstack([sp.kron(Dx, Idz), sp.kron(Idx, Dz)]).tocsr()
    A = M + c * (Div.T @ Div)
    perm = interleaved_permutation(Nx, Nz)
    A = A.tocsr()[perm][:, perm].tocsr()
    if return_parts:
        return A, M[perm][:, perm].tocsr(), Div[:, perm].tocsr()
    return A


def patch_indices(Nx, Nz, kind):
    if Nx < 3:
        raise ValueError('Nx >= 3 is needed')
    nzd = 2 * Nz
    nb = 8 * Nz - 2
    i = np.arange(Nx)[:, None]
    full = np.arange(nb)[None, :]
    vert = np.arange(nzd)[None, :]
    inter = np.arange(nzd, nb)[None, :]
    if kind == 'star':
        return np.hstack([((i - 1) % Nx) * nb + inter, i * nb + full])
    elif kind == 'column':
        return np.hstack([i * nb + full, ((i + 1) % Nx) * nb + vert])
    raise ValueError(kind)


class ASM:
    def __init__(self, A, Nx, Nz, kind):
        self.idx = patch_indices(Nx, Nz, kind)
        self.n = A.shape[0]
        A = A.tocsr()
        self.Ap = A[self.idx[0]][:, self.idx[0]].toarray()
        Apinv = sla.inv(self.Ap)
        self.Apinv = 0.5 * (Apinv + Apinv.T)
        self.flat = self.idx.ravel()

    def apply(self, r):
        y = r[self.idx] @ self.Apinv
        return np.bincount(self.flat, weights=y.ravel(), minlength=self.n)

    def dense(self):
        B = np.zeros((self.n, self.n))
        for i in range(self.idx.shape[0]):
            ii = self.idx[i]
            B[np.ix_(ii, ii)] += self.Apinv
        return B


def spectrum_direct(Nx, Nz, dx, dz, m, c, kind):
    A = assemble(Nx, Nz, dx, dz, m, c)
    B = ASM(A, Nx, Nz, kind).dense()
    Ad = A.toarray()
    L = np.linalg.cholesky(Ad)
    S = L.T @ B @ L
    return np.sort(np.linalg.eigvalsh(0.5 * (S + S.T)))


def spectrum_block_fourier(Nx, Nz, dx, dz, m, c, kind, thetas=None):
    '''Eigenvalues of B A by Fourier transform in x.  Returns (thetas, lam[ntheta, nb]).

    The blocks of A are taken from a small assembled matrix (4 columns, the same dx), so
    the cost does not depend on Nx.
    '''
    nzd = 2 * Nz
    nb = 8 * Nz - 2
    Ns = 4
    A = assemble(Ns, Nz, dx, dz, m, c).tocsr()
    rows = np.arange(nb, 2 * nb)
    Am = A[rows][:, np.arange(0, nb)].toarray()
    A0 = A[rows][:, rows].toarray()
    Ap = A[rows][:, np.arange(2 * nb, 3 * nb)].toarray()
    Binv = ASM(A, Ns, Nz, kind).Apinv
    npatch = Binv.shape[0]
    if thetas is None:
        thetas = 2.0 * np.pi * np.arange(Nx) / Nx
    lam = np.empty((len(thetas), nb))
    for q, th in enumerate(thetas):
        Ah = A0 + Ap * np.exp(1j * th) + Am * np.exp(-1j * th)
        Ah = 0.5 * (Ah + Ah.conj().T)
        E = np.zeros((npatch, nb), dtype=complex)
        if kind == 'star':       # [interior of block i-1, block i]
            E[np.arange(nb - nzd), np.arange(nzd, nb)] = np.exp(-1j * th)
            E[nb - nzd + np.arange(nb), np.arange(nb)] = 1.0
        else:                    # [block i, vertex part of block i+1]
            E[np.arange(nb), np.arange(nb)] = 1.0
            E[nb + np.arange(nzd), np.arange(nzd)] = np.exp(1j * th)
        Bh = E.conj().T @ Binv @ E
        L = np.linalg.cholesky(Ah)
        S = L.conj().T @ Bh @ L
        lam[q] = np.linalg.eigvalsh(0.5 * (S + S.conj().T))
    return np.asarray(thetas), lam


def spectral_summary(Nx, Nz, dx, dz, m, c, kind, frac=0.1):
    '''Same quantities as rt0_model.spectral_summary.  Only theta in [0, pi] is computed
    (the spectrum is symmetric under theta -> 2 pi - theta) and the multiplicities are
    accounted for.'''
    q = np.arange(Nx // 2 + 1)
    thetas = 2.0 * np.pi * q / Nx
    mult = np.where((q == 0) | ((Nx % 2 == 0) & (q == Nx // 2)), 1, 2)
    th, lam = spectrum_block_fourier(Nx, Nz, dx, dz, m, c, kind, thetas=thetas)
    w = np.repeat(mult[:, None], lam.shape[1], axis=1)
    lmax = lam.max()
    hf = thetas >= np.pi / 2 - 1e-12
    out = dict(
        ndof=int(w.sum()),
        lmin=lam.min(), lmax=lmax, cond=lmax / lam.min(),
        n_small=int(w[lam < frac * lmax].sum()),
        frac_small=float(w[lam < frac * lmax].sum() / w.sum()),
        ndof_hf=int(w[hf].sum()),
        lmin_hf=lam[hf].min(), cond_hf=lmax / lam[hf].min(),
        n_small_hf=int(w[hf][lam[hf] < frac * lmax].sum()),
        frac_small_hf=float(w[hf][lam[hf] < frac * lmax].sum() / w[hf].sum()),
    )
    # without the single smallest eigenvalue at theta = 0 if it is the uniform flow
    lam_nh = lam.copy()
    lam_nh[0, 0] = np.inf
    out['lmin_nh'] = lam_nh.min()
    out['cond_nh'] = lmax / lam_nh.min()
    return out


# ----------------------------------------------------------------------------
# Multigrid
# ----------------------------------------------------------------------------
def prolongation(Nxc, Nz):
    '''Natural embedding coarse -> fine (cells halved in x), interleaved numbering.'''
    Nxf = 2 * Nxc
    # CG2 periodic: coarse nodes 2I (vertex), 2I+1 (mid); fine nodes 4I, 4I+1, 4I+2, 4I+3
    rows, cols, vals = [], [], []
    nc = 2 * Nxc
    for I in range(Nxc):
        v0, mid, v1 = 2 * I, 2 * I + 1, (2 * I + 2) % nc
        rows += [4 * I]; cols += [v0]; vals += [1.0]
        rows += [4 * I + 2]; cols += [mid]; vals += [1.0]
        rows += [4 * I + 1] * 3; cols += [v0, mid, v1]; vals += [3.0 / 8.0, 3.0 / 4.0, -1.0 / 8.0]
        rows += [4 * I + 3] * 3; cols += [v0, mid, v1]; vals += [-1.0 / 8.0, 3.0 / 4.0, 3.0 / 8.0]
    Pcg = sp.csr_matrix((vals, (rows, cols)), shape=(2 * Nxf, nc))
    # DG1 orthonormal Legendre: coarse cell I -> fine cells 2I (left), 2I+1 (right)
    s2 = np.sqrt(2.0)
    s3 = np.sqrt(3.0)
    rows, cols, vals = [], [], []
    for I in range(Nxc):
        c0, c1 = 2 * I, 2 * I + 1
        fl0, fl1 = 2 * (2 * I), 2 * (2 * I) + 1
        fr0, fr1 = 2 * (2 * I + 1), 2 * (2 * I + 1) + 1
        rows += [fl0, fl0, fl1]; cols += [c0, c1, c1]; vals += [1.0 / s2, -s3 / (2.0 * s2), 1.0 / (2.0 * s2)]
        rows += [fr0, fr0, fr1]; cols += [c0, c1, c1]; vals += [1.0 / s2, s3 / (2.0 * s2), 1.0 / (2.0 * s2)]
    Pdg = sp.csr_matrix((vals, (rows, cols)), shape=(2 * Nxf, 2 * Nxc))
    Pfs = sp.block_diag([sp.kron(Pcg, sp.identity(2 * Nz)), sp.kron(Pdg, sp.identity(2 * Nz - 1))]).tocsr()
    pf = interleaved_permutation(Nxf, Nz)
    pc = interleaved_permutation(Nxc, Nz)
    return Pfs[pf][:, pc].tocsr()


class Multigrid:
    '''Same algorithm as rt0_model.Multigrid for the degree 2 element.'''

    def __init__(self, Nx, Nz, Lx, Lz, m, c, kind, nlevels=3, smooth_its=6, cycle='v'):
        assert Nx % (2 ** (nlevels - 1)) == 0
        self.nlevels = nlevels
        self.smooth_its = smooth_its
        self.cycle = cycle
        self.A, self.B, self.P = [], [], [None]
        dz = Lz / Nz
        for lev in range(nlevels):
            Nxl = Nx // (2 ** (nlevels - 1 - lev))
            A = assemble(Nxl, Nz, Lx / Nxl, dz, m, c)
            self.A.append(A)
            if lev == 0:
                self.B.append(None)
                self.coarse = spla.splu(A.tocsc())
            else:
                self.B.append(ASM(A, Nxl, Nz, kind))
                self.P.append(prolongation(Nxl // 2, Nz))
        self.n = self.A[-1].shape[0]

    def smooth(self, lev, b, x0):
        A = self.A[lev]
        return gmres_fixed(lambda v: A @ v, self.B[lev].apply, b, x0, self.smooth_its)

    def vcycle(self, lev, b, x0=None):
        if lev == 0:
            return self.coarse.solve(b)
        x = self.smooth(lev, b, x0)
        r = b - self.A[lev] @ x
        xc = self.vcycle(lev - 1, self.P[lev].T @ r)
        x = x + self.P[lev] @ xc
        return self.smooth(lev, b, x)

    def apply(self, b):
        L = self.nlevels - 1
        if self.cycle == 'v':
            return self.vcycle(L, b)
        bs = [None] * self.nlevels
        bs[L] = b
        for lev in range(L, 0, -1):
            bs[lev - 1] = self.P[lev].T @ bs[lev]
        x = self.vcycle(0, bs[0])
        for lev in range(1, self.nlevels):
            x = self.vcycle(lev, bs[lev], self.P[lev] @ x)
        return x


def mg_fgmres_iterations(Nx, Nz, Lx, Lz, m, c, kind, nlevels=3, smooth_its=6, cycle='v',
                         rtol=1e-8, seed=0, maxit=300, restart=30, return_history=False):
    mg = Multigrid(Nx, Nz, Lx, Lz, m, c, kind, nlevels=nlevels, smooth_its=smooth_its, cycle=cycle)
    A = mg.A[-1]
    rng = np.random.default_rng(seed)
    b = rng.standard_normal(mg.n)
    x, its, hist, conv = fgmres(lambda v: A @ v, mg.apply, b, rtol=rtol, maxit=maxit, restart=restart)
    true_res = np.linalg.norm(b - A @ x) / np.linalg.norm(b)
    if return_history:
        return its, conv, true_res, hist
    return its, conv, true_res
