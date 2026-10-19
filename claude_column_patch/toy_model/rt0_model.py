'''
Lowest-order toy model of the velocity block of the HDivSchurPC preconditioner.

Bilinear form (velocity block of field 1, linearised about rest, f = 0, b eliminated)

    a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w),
    m = 1 + (dt N)^2 / 4,     c = dt / (2 delta),

on the rectangle [0, L_x] x [0, L_z], periodic in x, u.n = 0 on top and bottom,
discretised with the lowest-order Raviart-Thomas element on a uniform mesh of
Nx x Nz rectangles of size dx x dz.

Degrees of freedom (point values of the normal component, NOT integrated fluxes):

    U[i, j]  u_x on the vertical facet above base vertex i in layer j,
             i = 0..Nx-1 (periodic), j = 0..Nz-1.      u_x is CG1 in x, DG0 in z.
    V[i, j]  u_z on the horizontal facet at height j*dz in cell column i,
             i = 0..Nx-1, j = 1..Nz-1 (top and bottom are zero). DG0 in x, CG1 in z.

Cell column i lies between base vertices i and i+1.

Global ("interleaved") numbering, block-circulant in x with block size nb = 2 Nz - 1:

    U[i, j] -> i*nb + j,            V[i, j] -> i*nb + Nz + (j - 1).

Patches (each patch is made of whole columns, as in Firedrake on an extruded mesh,
where the DMPlex is the base interval mesh):

    star   of base vertex i : V[i-1, :], U[i, :], V[i, :]          (3 Nz - 2 DoFs)
    column of base cell   i : U[i, :],   V[i, :], U[i+1, :]        (3 Nz - 1 DoFs)

Everything in this file is plain numpy / scipy.
'''
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.linalg as sla

PATCHES = ('star', 'column')


# ----------------------------------------------------------------------------
# Physical parameters
# ----------------------------------------------------------------------------
def physical_parameters(Lx, Lz, dt, K, N=1.0e-2):
    '''Return a dict with m, c, delta, C_M^2 and the regime, as in claude_paper_shift.py.

    delta = K * delta_thr, delta_thr = dt / (2 C_M^2), C_M^2 = max(L_x^2/(4 pi^2), m L_z^2/pi^2),
    hence c = dt / (2 delta) = C_M^2 / K.
    '''
    m = 1.0 + (dt * N) ** 2 / 4.0
    Ch2 = (Lx / (2.0 * np.pi)) ** 2
    Cv2 = (Lz / np.pi) ** 2
    CM2 = max(Ch2, m * Cv2)
    delta_thr = dt / (2.0 * CM2)
    delta = K * delta_thr
    c = dt / (2.0 * delta)
    eps = Lz / Lx
    return dict(m=m, c=c, delta=delta, delta_thr=delta_thr, CM2=CM2, eps=eps,
                regime='weak' if m * Cv2 <= Ch2 else 'strong',
                four_eps2_m=4.0 * eps ** 2 * m)


# ----------------------------------------------------------------------------
# 1D building blocks
# ----------------------------------------------------------------------------
def shift_matrix(n):
    '''Cyclic shift S with (S x)_i = x_{i+1 mod n}.'''
    rows = np.arange(n)
    cols = (rows + 1) % n
    return sp.csr_matrix((np.ones(n), (rows, cols)), shape=(n, n))


def mass_cg1_periodic(n, h):
    S = shift_matrix(n)
    return (h * (2.0 / 3.0 * sp.identity(n) + (S + S.T) / 6.0)).tocsr()


def diff_periodic(n):
    '''(Dx U)_i = U_{i+1} - U_i : vertex values -> cell column i (between vertices i, i+1).'''
    return (shift_matrix(n) - sp.identity(n)).tocsr()


def mass_cg1_dirichlet(nz, h):
    '''CG1 mass matrix in z restricted to the interior nodes j = 1..nz-1.'''
    n = nz - 1
    return (h * sp.diags([np.full(n - 1, 1.0 / 6.0), np.full(n, 2.0 / 3.0), np.full(n - 1, 1.0 / 6.0)],
                         [-1, 0, 1])).tocsr()


def diff_dirichlet(nz):
    '''(Dz V)_j = V_{j+1} - V_j for layer j = 0..nz-1 with V_0 = V_nz = 0; shape nz x (nz-1).'''
    rows, cols, vals = [], [], []
    for j in range(nz):
        if j + 1 <= nz - 1:          # V_{j+1}, interior index j
            rows.append(j); cols.append(j); vals.append(1.0)
        if j >= 1:                   # V_j, interior index j-1
            rows.append(j); cols.append(j - 1); vals.append(-1.0)
    return sp.csr_matrix((vals, (rows, cols)), shape=(nz, nz - 1))


# ----------------------------------------------------------------------------
# Numbering
# ----------------------------------------------------------------------------
def fieldsplit_to_interleaved(Nx, Nz):
    '''perm[g] = index in the field-split ordering [all U (i*Nz+j); all V (i*(Nz-1)+j-1)]
    of the DoF with interleaved global index g.'''
    nb = 2 * Nz - 1
    nU = Nx * Nz
    perm = np.empty(Nx * nb, dtype=np.int64)
    for i in range(Nx):
        perm[i * nb: i * nb + Nz] = i * Nz + np.arange(Nz)
        perm[i * nb + Nz: (i + 1) * nb] = nU + i * (Nz - 1) + np.arange(Nz - 1)
    return perm


def permute(A, perm):
    A = A.tocsr()
    return A[perm][:, perm].tocsr()


# ----------------------------------------------------------------------------
# Assembly
# ----------------------------------------------------------------------------
def assemble(Nx, Nz, dx, dz, m, c, return_parts=False):
    '''Assemble A = M + c * dx dz * D^T D in the interleaved numbering.

    M : consistent RT0 mass matrix weighted with (1, m), D : cell-wise divergence.
    If return_parts: return (A, M_unweighted_by_c, Dmat) with Dmat mapping to cells (i*Nz + j).
    '''
    Mx = mass_cg1_periodic(Nx, dx)
    Dx = diff_periodic(Nx)
    Mz = mass_cg1_dirichlet(Nz, dz)
    Dz = diff_dirichlet(Nz)
    Ix = sp.identity(Nx)
    Izc = sp.identity(Nz)
    M_UU = sp.kron(Mx, dz * Izc)
    M_VV = sp.kron(dx * Ix, m * Mz)
    Mfs = sp.block_diag([M_UU, M_VV]).tocsr()
    Dfs = sp.hstack([sp.kron(Dx, Izc) / dx, sp.kron(Ix, Dz) / dz]).tocsr()
    Afs = Mfs + c * dx * dz * (Dfs.T @ Dfs)
    perm = fieldsplit_to_interleaved(Nx, Nz)
    A = permute(Afs, perm)
    if return_parts:
        M = permute(Mfs, perm)
        D = Dfs[:, perm].tocsr()
        return A, M, D
    return A


def assemble_by_quadrature(Nx, Nz, dx, dz, m, c, nq=3):
    '''Independent element-by-element assembly with Gauss quadrature (validation only).

    Local RT0 basis on the reference rectangle [0,dx] x [0,dz]:
        left facet  : ((1 - x/dx), 0)      right facet : (x/dx, 0)
        bottom facet: (0, (1 - z/dz))      top facet   : (0, z/dz)
    '''
    nb = 2 * Nz - 1
    n = Nx * nb
    A = np.zeros((n, n))
    xg, wg = np.polynomial.legendre.leggauss(nq)
    xg = 0.5 * (xg + 1.0)
    wg = 0.5 * wg

    def gU(i, j):
        return (i % Nx) * nb + j

    def gV(i, j):  # j = 0..Nz, boundary -> -1
        if j == 0 or j == Nz:
            return -1
        return (i % Nx) * nb + Nz + (j - 1)

    for i in range(Nx):
        for j in range(Nz):
            dofs = [gU(i, j), gU(i + 1, j), gV(i, j), gV(i, j + 1)]
            Ae = np.zeros((4, 4))
            for a, wa in zip(xg, wg):
                for b, wb in zip(xg, wg):
                    w = wa * wb * dx * dz
                    # basis values (component x, component z) and divergence
                    phi = np.array([[1 - a, 0.0], [a, 0.0], [0.0, 1 - b], [0.0, b]])
                    div = np.array([-1.0 / dx, 1.0 / dx, -1.0 / dz, 1.0 / dz])
                    Ae += w * (np.outer(phi[:, 0], phi[:, 0]) + m * np.outer(phi[:, 1], phi[:, 1])
                               + c * np.outer(div, div))
            for p in range(4):
                if dofs[p] < 0:
                    continue
                for q in range(4):
                    if dofs[q] < 0:
                        continue
                    A[dofs[p], dofs[q]] += Ae[p, q]
    return A


# ----------------------------------------------------------------------------
# Patches
# ----------------------------------------------------------------------------
def patch_indices(Nx, Nz, kind):
    '''Integer array idx[i, :] with the global DoFs of patch i (i = 0..Nx-1).

    The local ordering is the same for every i, so that the patch matrix is the same
    for all patches (translation invariance).
    star   : [V[i-1, :], U[i, :], V[i, :]]
    column : [U[i, :],   V[i, :], U[i+1, :]]
    '''
    nb = 2 * Nz - 1
    i = np.arange(Nx)[:, None]
    jU = np.arange(Nz)[None, :]
    jV = Nz + np.arange(Nz - 1)[None, :]
    if kind == 'star':
        idx = np.hstack([((i - 1) % Nx) * nb + jV, i * nb + jU, i * nb + jV])
    elif kind == 'column':
        idx = np.hstack([i * nb + jU, i * nb + jV, ((i + 1) % Nx) * nb + jU])
    else:
        raise ValueError(kind)
    if Nx < 3:
        raise ValueError('Nx >= 3 is needed for well defined patches')
    return idx


class ASM:
    '''Additive Schwarz  B = sum_i R_i^T A_i^{-1} R_i  with exact patch solves.

    All patch matrices are equal (uniform periodic mesh), so a single dense inverse is
    stored and applied to all patches at once.
    '''

    def __init__(self, A, Nx, Nz, kind, check_translation_invariance=False):
        self.kind = kind
        self.idx = patch_indices(Nx, Nz, kind)
        self.n = A.shape[0]
        A = A.tocsr()
        Ap = A[self.idx[0]][:, self.idx[0]].toarray()
        if check_translation_invariance:
            for i in range(1, Nx):
                Ai = A[self.idx[i]][:, self.idx[i]].toarray()
                assert np.allclose(Ai, Ap, rtol=1e-13, atol=0.0)
        self.Ap = Ap
        self.Apinv = sla.inv(Ap)
        self.Apinv = 0.5 * (self.Apinv + self.Apinv.T)
        self.flat = self.idx.ravel()

    def apply(self, r):
        y = r[self.idx] @ self.Apinv          # (Nx, np); Apinv symmetric
        return np.bincount(self.flat, weights=y.ravel(), minlength=self.n)

    def dense(self):
        B = np.zeros((self.n, self.n))
        for i in range(self.idx.shape[0]):
            ii = self.idx[i]
            B[np.ix_(ii, ii)] += self.Apinv
        return B


# ----------------------------------------------------------------------------
# Spectra: direct (dense), block Fourier in x, and Fourier in x + sine/cosine in z
# ----------------------------------------------------------------------------
def spectrum_direct(Nx, Nz, dx, dz, m, c, kind):
    '''All eigenvalues of B A from dense matrices (generalised symmetric eigenproblem).'''
    A = assemble(Nx, Nz, dx, dz, m, c).toarray()
    B = ASM(sp.csr_matrix(A), Nx, Nz, kind).dense()
    # eigenvalues of B A = eigenvalues of L^T B L with A = L L^T
    L = np.linalg.cholesky(A)
    S = L.T @ B @ L
    S = 0.5 * (S + S.T)
    return np.sort(np.linalg.eigvalsh(S))


def _block_symbols(Nx, Nz, dx, dz, m, c):
    '''Blocks A_{-1}, A_0, A_{+1} (size nb x nb) of the block-circulant matrix A, taken
    from the assembled matrix.'''
    nb = 2 * Nz - 1
    A = assemble(Nx, Nz, dx, dz, m, c).tocsr()
    i0 = 1  # a generic block row
    rows = np.arange(i0 * nb, (i0 + 1) * nb)
    Am = A[rows][:, np.arange((i0 - 1) * nb, i0 * nb)].toarray()
    A0 = A[rows][:, rows].toarray()
    Ap = A[rows][:, np.arange((i0 + 1) * nb, (i0 + 2) * nb)].toarray()
    return A, Am, A0, Ap


def spectrum_block_fourier(Nx, Nz, dx, dz, m, c, kind, return_theta=False):
    '''Eigenvalues of B A by Fourier transform in x (blocks of size 2 Nz - 1).

    For theta = 2 pi q / Nx the mode is x_i = xhat * exp(1j theta i).
    '''
    assert Nx >= 4
    nb = 2 * Nz - 1
    A, Am, A0, Ap = _block_symbols(Nx, Nz, dx, dz, m, c)
    asm = ASM(A, Nx, Nz, kind)
    Binv_loc = asm.Apinv
    npatch = Binv_loc.shape[0]
    thetas = 2.0 * np.pi * np.arange(Nx) / Nx
    lam = np.empty((Nx, nb))
    for q, th in enumerate(thetas):
        Ah = A0 + Ap * np.exp(1j * th) + Am * np.exp(-1j * th)
        Ah = 0.5 * (Ah + Ah.conj().T)
        # extraction E(theta): block vector (U, V) -> patch vector
        E = np.zeros((npatch, nb), dtype=complex)
        if kind == 'star':       # [V[i-1], U[i], V[i]]
            E[0:Nz - 1, Nz:] = np.exp(-1j * th) * np.eye(Nz - 1)
            E[Nz - 1:2 * Nz - 1, :Nz] = np.eye(Nz)
            E[2 * Nz - 1:, Nz:] = np.eye(Nz - 1)
        else:                    # [U[i], V[i], U[i+1]]
            E[0:Nz, :Nz] = np.eye(Nz)
            E[Nz:2 * Nz - 1, Nz:] = np.eye(Nz - 1)
            E[2 * Nz - 1:, :Nz] = np.exp(1j * th) * np.eye(Nz)
        Bh = E.conj().T @ Binv_loc @ E
        Bh = 0.5 * (Bh + Bh.conj().T)
        L = np.linalg.cholesky(Ah)
        S = L.conj().T @ Bh @ L
        S = 0.5 * (S + S.conj().T)
        lam[q] = np.linalg.eigvalsh(S)
    if return_theta:
        return thetas, lam
    return np.sort(lam.ravel())


def symbols_2x2(Nx, Nz, dx, dz, m, c, kind, thetas=None):
    '''Eigenvalues of B A from the fully diagonalised symbols.

    Fourier in x (theta) and the singular vectors of Dz in z (index l = 0..Nz-1):
        Dz v_l = sigma_l u_l,  v_l(j) = sin(l pi j / Nz),  u_l(j) = cos(l pi (j + 1/2) / Nz),
        sigma_l = 2 sin(l pi / (2 Nz)),  mu_l = sigma_l^2.
    For l >= 1 the symbol of A is 2 x 2 (one U and one V amplitude), for l = 0 it is a
    scalar (U only).

    Returns (thetas, lam1, lam0): lam1 has shape (ntheta, Nz-1, 2) (sorted ascending in the
    last index), lam0 has shape (ntheta,).
    '''
    if thetas is None:
        thetas = 2.0 * np.pi * np.arange(Nx) / Nx
    thetas = np.asarray(thetas, dtype=float)
    th = thetas[:, None]
    l = np.arange(1, Nz)[None, :]
    sig = 2.0 * np.sin(l * np.pi / (2.0 * Nz))
    mu = sig ** 2
    e = np.exp(1j * th)
    one = np.ones_like(th * sig)

    aUU = dx * dz * (2.0 + np.cos(th)) / 3.0 + c * dz / dx * (2.0 - 2.0 * np.cos(th))
    aVV = m * dx * dz * (1.0 - mu / 6.0) + c * dx / dz * mu
    aUV = c * np.conj(e - 1.0) * sig
    Ah = np.empty(one.shape + (2, 2), dtype=complex)
    Ah[..., 0, 0] = aUU * one
    Ah[..., 0, 1] = aUV
    Ah[..., 1, 0] = np.conj(aUV)
    Ah[..., 1, 1] = aVV * one

    dU = 2.0 * dx * dz / 3.0 + 2.0 * c * dz / dx      # diagonal U entry of a patch
    oU = dx * dz / 6.0 - c * dz / dx                  # coupling of the two U columns of a patch
    cs = c * sig * one
    Ap = np.zeros(one.shape + (3, 3))
    E = np.zeros(one.shape + (3, 2), dtype=complex)
    if kind == 'column':        # (UL, V, UR)
        Ap[..., 0, 0] = dU; Ap[..., 2, 2] = dU
        Ap[..., 0, 2] = oU; Ap[..., 2, 0] = oU
        Ap[..., 1, 1] = aVV * one
        Ap[..., 0, 1] = -cs; Ap[..., 1, 0] = -cs
        Ap[..., 2, 1] = cs; Ap[..., 1, 2] = cs
        E[..., 0, 0] = 1.0
        E[..., 1, 1] = 1.0
        E[..., 2, 0] = e * one
    elif kind == 'star':        # (VL, U, VR)
        Ap[..., 0, 0] = aVV * one; Ap[..., 2, 2] = aVV * one
        Ap[..., 1, 1] = dU
        Ap[..., 0, 1] = cs; Ap[..., 1, 0] = cs
        Ap[..., 1, 2] = -cs; Ap[..., 2, 1] = -cs
        E[..., 0, 1] = np.conj(e) * one
        E[..., 1, 0] = 1.0
        E[..., 2, 1] = 1.0
    else:
        raise ValueError(kind)
    Apinv = np.linalg.inv(Ap)
    EH = np.conj(np.swapaxes(E, -1, -2))
    Bh = EH @ Apinv @ E
    L = np.linalg.cholesky(Ah)
    LH = np.conj(np.swapaxes(L, -1, -2))
    S = LH @ Bh @ L
    S = 0.5 * (S + np.conj(np.swapaxes(S, -1, -2)))
    lam1 = np.linalg.eigvalsh(S)

    # l = 0 : U only
    th0 = thetas
    aUU0 = dx * dz * (2.0 + np.cos(th0)) / 3.0 + c * dz / dx * (2.0 - 2.0 * np.cos(th0))
    if kind == 'column':
        Ap0 = np.array([[dU, oU], [oU, dU]])
        Ap0inv = np.linalg.inv(Ap0)
        e0 = np.exp(1j * th0)
        # E = (1, e)^T
        b0 = (Ap0inv[0, 0] + Ap0inv[1, 1] + 2.0 * Ap0inv[0, 1] * np.cos(th0))
        lam0 = aUU0 * b0
    else:
        lam0 = aUU0 / dU
    return thetas, lam1, lam0


def spectrum_2x2(Nx, Nz, dx, dz, m, c, kind):
    th, lam1, lam0 = symbols_2x2(Nx, Nz, dx, dz, m, c, kind)
    return np.sort(np.concatenate([lam1.ravel(), lam0.ravel()]))


def spectral_summary(Nx, Nz, dx, dz, m, c, kind, frac=0.1):
    '''Summary of the spectrum of B A from the 2 x 2 symbols.

    all modes            : lmin, lmax, cond, number of eigenvalues below frac * lmax
    no uniform flow      : lmin_nh, cond_nh = lmax / lmin_nh, all modes except (theta = 0, l = 0)
    high frequency in x  : the same restricted to pi/2 <= theta <= 3 pi/2, i.e. the modes
                           that the next coarser level (coarsening in x by 2) cannot represent
    '''
    th, lam1, lam0 = symbols_2x2(Nx, Nz, dx, dz, m, c, kind)
    lam_all = np.concatenate([lam1.ravel(), lam0.ravel()])
    hf = (th >= np.pi / 2 - 1e-12) & (th <= 3 * np.pi / 2 + 1e-12)
    lam_hf = np.concatenate([lam1[hf].ravel(), lam0[hf].ravel()])
    lmax = lam_all.max()
    # all modes except the uniform horizontal flow u = (1, 0) (theta = 0, l = 0), which is the
    # one divergence-free function that is not a sum of local divergence-free functions
    lam_nh = np.concatenate([lam1.ravel(), lam0[1:]])
    assert th[0] == 0.0
    out = dict(
        ndof=lam_all.size,
        lmin=lam_all.min(), lmax=lmax, cond=lmax / lam_all.min(),
        n_small=int(np.sum(lam_all < frac * lmax)),
        frac_small=float(np.mean(lam_all < frac * lmax)),
        lam_uniform_flow=lam0[0],
        lmin_nh=lam_nh.min(), cond_nh=lmax / lam_nh.min(),
        ndof_hf=lam_hf.size,
        lmin_hf=lam_hf.min(), lmax_hf=lam_hf.max(), cond_hf=lmax / lam_hf.min(),
        n_small_hf=int(np.sum(lam_hf < frac * lmax)),
        frac_small_hf=float(np.mean(lam_hf < frac * lmax)),
    )
    # where is the smallest eigenvalue attained?
    k = np.unravel_index(np.argmin(lam1[..., 0]), lam1[..., 0].shape)
    if lam1[..., 0].min() <= lam0.min():
        out['argmin_theta'] = th[k[0]]
        out['argmin_l'] = int(k[1] + 1)
    else:
        out['argmin_theta'] = th[np.argmin(lam0)]
        out['argmin_l'] = 0
    return out


# ----------------------------------------------------------------------------
# Multigrid (coarsening in x only)
# ----------------------------------------------------------------------------
def prolongation(Nxc, Nz):
    '''Natural embedding of the coarse RT0 space (Nxc columns) into the fine one (2 Nxc
    columns), interleaved numbering on both levels.

    U (CG1 in x): fine vertex 2i = coarse vertex i, fine vertex 2i+1 = mean of coarse i, i+1.
    V (DG0 in x): fine columns 2i and 2i+1 = coarse column i.
    '''
    Nxf = 2 * Nxc
    rows, cols, vals = [], [], []
    for i in range(Nxc):
        rows += [2 * i, 2 * i + 1, 2 * i + 1]
        cols += [i, i, (i + 1) % Nxc]
        vals += [1.0, 0.5, 0.5]
    Pcg = sp.csr_matrix((vals, (rows, cols)), shape=(Nxf, Nxc))
    rows, cols, vals = [], [], []
    for i in range(Nxc):
        rows += [2 * i, 2 * i + 1]
        cols += [i, i]
        vals += [1.0, 1.0]
    Pdg = sp.csr_matrix((vals, (rows, cols)), shape=(Nxf, Nxc))
    Pfs = sp.block_diag([sp.kron(Pcg, sp.identity(Nz)), sp.kron(Pdg, sp.identity(Nz - 1))]).tocsr()
    pf = fieldsplit_to_interleaved(Nxf, Nz)
    pc = fieldsplit_to_interleaved(Nxc, Nz)
    return Pfs[pf][:, pc].tocsr()


def gmres_fixed(Aop, Bop, b, x0, nits):
    '''nits iterations of left-preconditioned GMRES (no convergence test), as the PETSc
    smoother  -mg_levels_ksp_type gmres -mg_levels_ksp_max_it nits  (norm type none).
    Minimises || B (b - A x) ||_2 over x0 + K_nits(B A, B r0).'''
    if x0 is None:
        x0 = np.zeros_like(b)
        r = Bop(b)
    else:
        r = Bop(b - Aop(x0))
    beta = np.linalg.norm(r)
    if beta == 0.0:
        return x0
    n = b.size
    Vk = np.empty((nits + 1, n))
    H = np.zeros((nits + 1, nits))
    Vk[0] = r / beta
    k_used = nits
    for k in range(nits):
        w = Bop(Aop(Vk[k]))
        for j in range(k + 1):
            H[j, k] = Vk[j] @ w
            w = w - H[j, k] * Vk[j]
        # one re-orthogonalisation pass for safety
        for j in range(k + 1):
            hc = Vk[j] @ w
            H[j, k] += hc
            w = w - hc * Vk[j]
        H[k + 1, k] = np.linalg.norm(w)
        if H[k + 1, k] <= 1e-14 * beta:
            k_used = k + 1
            break
        Vk[k + 1] = w / H[k + 1, k]
    rhs = np.zeros(k_used + 1)
    rhs[0] = beta
    y, *_ = np.linalg.lstsq(H[:k_used + 1, :k_used], rhs, rcond=None)
    return x0 + y @ Vk[:k_used]


class Multigrid:
    '''Geometric multigrid in x for A, levels 0 (coarsest) .. L (finest).

    smoother   : smooth_its iterations of GMRES preconditioned by the additive Schwarz
                 method with the given patch (pre and post smoothing, same number)
    coarse     : sparse LU
    cycle      : 'v'    one V cycle
                 'full' PETSc -pc_mg_type full (restrict the right-hand side to all levels,
                        solve on the coarsest level, then on every finer level interpolate
                        and do one V cycle)
    operators  : rediscretisation on every level (identical to Galerkin here, see validation)
    '''

    def __init__(self, Nx, Nz, Lx, Lz, m, c, kind, nlevels=3, smooth_its=6, cycle='v',
                 smoother='gmres', omega=None):
        assert Nx % (2 ** (nlevels - 1)) == 0
        self.nlevels = nlevels
        self.kind = kind
        self.smooth_its = smooth_its
        self.cycle = cycle
        self.smoother = smoother
        self.omega = omega
        self.Nz = Nz
        self.A, self.B, self.P, self.Nxs = [], [], [None], []
        dz = Lz / Nz
        for lev in range(nlevels):
            Nxl = Nx // (2 ** (nlevels - 1 - lev))
            self.Nxs.append(Nxl)
            A = assemble(Nxl, Nz, Lx / Nxl, dz, m, c)
            self.A.append(A)
            if lev == 0:
                self.B.append(None)
                self.coarse = spla.splu(A.tocsc())
            else:
                self.B.append(ASM(A, Nxl, Nz, kind))
                self.P.append(prolongation(self.Nxs[lev - 1], Nz))
        self.n = self.A[-1].shape[0]
        self.n_smooth_applications = 0

    def smooth(self, lev, b, x0):
        A = self.A[lev]
        B = self.B[lev]
        if self.smoother == 'gmres':
            self.n_smooth_applications += self.smooth_its
            return gmres_fixed(lambda v: A @ v, B.apply, b, x0, self.smooth_its)
        elif self.smoother == 'richardson':
            x = np.zeros_like(b) if x0 is None else x0.copy()
            for _ in range(self.smooth_its):
                x = x + self.omega * B.apply(b - A @ x)
            return x
        raise ValueError(self.smoother)

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
        elif self.cycle == 'full':
            bs = [None] * self.nlevels
            bs[L] = b
            for lev in range(L, 0, -1):
                bs[lev - 1] = self.P[lev].T @ bs[lev]
            x = self.vcycle(0, bs[0])
            for lev in range(1, self.nlevels):
                x = self.vcycle(lev, bs[lev], self.P[lev] @ x)
            return x
        raise ValueError(self.cycle)


def fgmres(Aop, Mop, b, rtol=1e-8, atol=0.0, restart=30, maxit=300):
    '''Right-preconditioned flexible GMRES, zero initial guess.

    Returns (x, iterations, history of true-recurrence residual norms, converged).
    The iteration count is the number of preconditioner applications.
    '''
    n = b.size
    x = np.zeros(n)
    bnorm = np.linalg.norm(b)
    tol = max(rtol * bnorm, atol)
    hist = [bnorm]
    its = 0
    r = b.copy()
    beta = bnorm
    if beta <= tol:
        return x, 0, hist, True
    while its < maxit:
        V = np.empty((restart + 1, n))
        Z = np.empty((restart, n))
        H = np.zeros((restart + 1, restart))
        V[0] = r / beta
        g = np.zeros(restart + 1)
        g[0] = beta
        cs = np.zeros(restart)
        sn = np.zeros(restart)
        k_used = 0
        converged = False
        for k in range(restart):
            Z[k] = Mop(V[k])
            w = Aop(Z[k])
            for j in range(k + 1):
                H[j, k] = V[j] @ w
                w = w - H[j, k] * V[j]
            for j in range(k + 1):
                hc = V[j] @ w
                H[j, k] += hc
                w = w - hc * V[j]
            H[k + 1, k] = np.linalg.norm(w)
            if H[k + 1, k] > 0.0:
                V[k + 1] = w / H[k + 1, k]
            # Givens rotations
            for j in range(k):
                t = cs[j] * H[j, k] + sn[j] * H[j + 1, k]
                H[j + 1, k] = -sn[j] * H[j, k] + cs[j] * H[j + 1, k]
                H[j, k] = t
            den = np.hypot(H[k, k], H[k + 1, k])
            cs[k] = H[k, k] / den
            sn[k] = H[k + 1, k] / den
            H[k, k] = den
            H[k + 1, k] = 0.0
            g[k + 1] = -sn[k] * g[k]
            g[k] = cs[k] * g[k]
            its += 1
            k_used = k + 1
            res = abs(g[k + 1])
            hist.append(res)
            if res <= tol or its >= maxit:
                converged = res <= tol
                break
        y = sla.solve_triangular(H[:k_used, :k_used], g[:k_used])
        x = x + y @ Z[:k_used]
        r = b - Aop(x)
        beta = np.linalg.norm(r)
        if converged or beta <= tol:
            return x, its, hist, True
    return x, its, hist, False


def mg_fgmres_iterations(Nx, Nz, Lx, Lz, m, c, kind, nlevels=3, smooth_its=6, cycle='v',
                         rtol=1e-8, seed=0, maxit=300, restart=30, smoother='gmres', omega=None,
                         return_history=False):
    '''FGMRES iterations for A x = b, b random, preconditioned by one multigrid cycle.'''
    mg = Multigrid(Nx, Nz, Lx, Lz, m, c, kind, nlevels=nlevels, smooth_its=smooth_its,
                   cycle=cycle, smoother=smoother, omega=omega)
    A = mg.A[-1]
    rng = np.random.default_rng(seed)
    b = rng.standard_normal(mg.n)
    x, its, hist, conv = fgmres(lambda v: A @ v, mg.apply, b, rtol=rtol, maxit=maxit, restart=restart)
    true_res = np.linalg.norm(b - A @ x) / np.linalg.norm(b)
    if return_history:
        return its, conv, true_res, hist
    return its, conv, true_res
