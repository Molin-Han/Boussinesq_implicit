'''
Exact modal reduction of the velocity operator of field 1 (linearised about rest, f = 0,
buoyancy eliminated)

    a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w)

on the extruded mesh of the solver (periodic in x, N_x columns of width dx, N_z layers,
u.n = 0 on top and bottom, RT_k on rectangles: u_x in CG_k(x) x DG_{k-1}(z),
u_z in DG_{k-1}(x) x CG_k(z)).

Because every patch of the solver (star or vanka) is a union of WHOLE columns, and the
multigrid hierarchy coarsens in x only, the operator, both additive Schwarz methods and
the transfer operators are block diagonal in the basis of vertical modes

    K_z chi_n = lambda_n^2 M_z chi_n ,   chi_n in CG_k(z), chi_n = 0 at z = 0, L_z,
    psi_n = chi_n' / lambda_n  in DG_{k-1}(z),   psi_0 = const  (the barotropic mode),

and for the mode n the problem is the ONE-dimensional problem in x

    a_n((U, W), (U~, W~)) = (U, U~) + m (W, W~) + c (U' + lambda_n W, U~' + lambda_n W~),
    U in CG_k(x),  W in DG_{k-1}(x)            (for n = 0 there is no W).

After scaling W by sqrt(m) and lengths by dx this problem has exactly two parameters,

    alpha = c / dx^2                (horizontal grad-div against the u mass),
    b_n   = c lambda_n^2 / m        (vertical grad-div against the w mass),

and eliminating W cell by cell gives the 1D reaction-diffusion operator
    S_n = M + c_eff K,   c_eff = c / (1 + b_n).

This module builds the 1D matrices (k = 1 or 2), the two patch decompositions
    star  of the base vertex i   : {U vertex i, U interior and W of cells i-1/2 and i+1/2}
    vanka (closure of base cell) : {U vertices i and i+1, U interior and W of cell i+1/2}
the x-coarsening prolongation, and the tools used by the other scripts in this folder.
'''
import numpy as np
import scipy.linalg as sla


# ----------------------------------------------------------------------------------
# reference matrices on the unit cell
# ----------------------------------------------------------------------------------
def reference(k):
    '''
    Returns (Mu, Ku, C, nw) on the unit cell.
      u local dofs: k = 1: (left vertex, right vertex)
                    k = 2: (left vertex, midpoint, right vertex)   (nodal Lagrange basis)
      w local dofs: orthonormal Legendre basis of DG_{k-1}, so that the w mass is the identity.
      C[a, j] = int psi_a phi_j' dxi.
    '''
    if k == 1:
        Mu = np.array([[2., 1.], [1., 2.]]) / 6.
        Ku = np.array([[1., -1.], [-1., 1.]])
        C = np.array([[-1., 1.]])
    elif k == 2:
        Mu = np.array([[4., 2., -1.], [2., 16., 2.], [-1., 2., 4.]]) / 30.
        Ku = np.array([[7., -8., 1.], [-8., 16., -8.], [1., -8., 7.]]) / 3.
        s3 = np.sqrt(3.)
        C = np.array([[-1., 0., 1.],
                      [2. * s3 / 3., -4. * s3 / 3., 2. * s3 / 3.]])
    else:
        raise ValueError('k must be 1 or 2')
    assert np.allclose(C.T @ C, Ku)
    return Mu, Ku, C, C.shape[0]


class Level:
    '''
    One level of the 1D problem: N periodic cells of width dx.
    Global numbering:
        k = 1:  u vertex i -> i                       (N)
                w cell i   -> N + i                   (N)
        k = 2:  u vertex i -> i                       (N)
                u mid  i   -> N + i                   (N)     cell i = (vertex i, vertex i+1)
                w cell i   -> 2N + 2 i + (0, 1)       (2N)
    '''

    def __init__(self, N, dx, k):
        self.N, self.dx, self.k = N, dx, k
        Mu, Ku, C, nw = reference(k)
        self.nw = nw
        self.nu = N * k
        self.n = self.nu + N * nw
        nu = self.nu
        self.Mu = np.zeros((nu, nu))
        self.Ku = np.zeros((nu, nu))
        self.C = np.zeros((N * nw, nu))
        for i in range(N):
            ud = self.cell_u_dofs(i)
            wd = self.cell_w_dofs(i, local=True)
            self.Mu[np.ix_(ud, ud)] += dx * Mu
            self.Ku[np.ix_(ud, ud)] += Ku / dx
            self.C[np.ix_(wd, ud)] += C
        self.Mw = dx * np.eye(N * nw)

    def cell_u_dofs(self, i):
        N = self.N
        if self.k == 1:
            return [i, (i + 1) % N]
        return [i, N + i, (i + 1) % N]

    def cell_w_dofs(self, i, local=False):
        off = 0 if local else self.nu
        return [off + self.nw * i + a for a in range(self.nw)]

    def cell_interior_dofs(self, i):
        '''u interior (k = 2) and w dofs of cell i, global numbering.'''
        d = [] if self.k == 1 else [self.N + i]
        return d + self.cell_w_dofs(i)

    # ------------------------------------------------------------------
    def matrix(self, c, lam, m=1.0, with_w=True):
        '''
        a_n in physical units; lam = lambda_n; unknowns (U, W), W NOT rescaled.
        with_w = False returns the barotropic (n = 0) operator M + c K on U only.
        '''
        Auu = self.Mu + c * self.Ku
        if not with_w:
            return Auu
        Awu = c * lam * self.C
        Aww = (m + c * lam ** 2) * self.Mw
        return np.block([[Auu, Awu.T], [Awu, Aww]])

    def schur_u(self, c, lam, m=1.0):
        '''Exact Schur complement on U: M + c_eff K with c_eff = c / (1 + c lam^2 / m).'''
        return self.Mu + c / (1.0 + c * lam ** 2 / m) * self.Ku

    # ------------------------------------------------------------------
    def patches(self, kind, with_w=True):
        N = self.N
        out = []
        for i in range(N):
            if kind == 'star':
                # star of vertex i: the vertex dof and the interiors of cells i-1 and i
                d = [i]
                for cell in ((i - 1) % N, i):
                    if self.k == 2:
                        d.append(N + cell)
                    if with_w:
                        d += self.cell_w_dofs(cell)
            elif kind == 'vanka':
                # closure of cell i: both vertices and the interior of cell i
                d = [i, (i + 1) % N]
                if self.k == 2:
                    d.append(N + i)
                if with_w:
                    d += self.cell_w_dofs(i)
            else:
                raise ValueError(kind)
            out.append(np.array(sorted(set(d))))
        return out

    def ndof(self, with_w=True):
        return self.n if with_w else self.nu


def asm_matrix(A, patches):
    '''B = sum_p R_p^T A_p^{-1} R_p as a dense matrix.'''
    B = np.zeros_like(A)
    for p in patches:
        B[np.ix_(p, p)] += np.linalg.inv(A[np.ix_(p, p)])
    return B


def prolongation(coarse, fine, with_w=True):
    '''
    Natural embedding of the coarse space (cells of width 2 dx) into the fine one.
    Coarse cell I = fine cells 2I and 2I+1; coarse vertex I = fine vertex 2I.
    '''
    k = fine.k
    Nc, Nf = coarse.N, fine.N
    assert Nf == 2 * Nc
    P = np.zeros((fine.ndof(with_w), coarse.ndof(with_w)))
    # --- u: evaluate coarse basis functions at the fine nodes
    if k == 1:
        for I in range(Nc):
            P[2 * I, I] += 1.0
            P[2 * I + 1, I] += 0.5
            P[2 * I + 1, (I + 1) % Nc] += 0.5
    else:
        def phi(xi):
            return np.array([(1 - xi) * (1 - 2 * xi), 4 * xi * (1 - xi), xi * (2 * xi - 1)])
        for I in range(Nc):
            cd = [I, Nc + I, (I + 1) % Nc]
            # fine nodes inside coarse cell I (excluding the right vertex, owned by the next cell)
            # fine vertex 2I at xi = 0, fine mid of cell 2I at 1/4, fine vertex 2I+1 at 1/2,
            # fine mid of cell 2I+1 at 3/4
            for fd, xi in ((2 * I, 0.0), (Nf + 2 * I, 0.25), (2 * I + 1, 0.5), (Nf + 2 * I + 1, 0.75)):
                P[fd, cd] += phi(xi)
    if with_w:
        # --- w: coarse orthonormal Legendre functions expressed in the fine ones.
        # coarse psi_0 = 1, psi_1 = sqrt(3)(2 X - 1), X in (0,1) coarse coordinate.
        # On the left fine cell X = xi/2, on the right X = (1 + xi)/2.
        for I in range(Nc):
            for half, fc in ((0, 2 * I), (1, 2 * I + 1)):
                fw = fine.cell_w_dofs(fc)
                cw = coarse.cell_w_dofs(I)
                if k == 1:
                    P[fw[0], cw[0]] = 1.0
                else:
                    # psi_0^c = psi_0^f
                    P[fw[0], cw[0]] = 1.0
                    # psi_1^c = sqrt3 (2X - 1); left: sqrt3 (xi - 1) = 0.5 sqrt3 (2 xi - 1) - 0.5 sqrt3
                    #                            right: sqrt3 xi    = 0.5 sqrt3 (2 xi - 1) + 0.5 sqrt3
                    P[fw[1], cw[1]] = 0.5
                    P[fw[0], cw[1]] = -0.5 * np.sqrt(3.) if half == 0 else 0.5 * np.sqrt(3.)
    return P


# ----------------------------------------------------------------------------------
# vertical eigenvalues lambda_n^2 of CG_k with homogeneous Dirichlet conditions
# ----------------------------------------------------------------------------------
def vertical_eigenvalues(Nz, Lz, k):
    '''
    Generalised eigenvalues lambda_n^2 of K_z chi = lambda^2 M_z chi for CG_k on Nz
    uniform layers of (0, Lz) with chi(0) = chi(Lz) = 0.  Returned sorted, length k Nz - 1.
    '''
    dz = Lz / Nz
    Mu, Ku, _, _ = reference(k)
    n = k * Nz + 1
    M = np.zeros((n, n))
    K = np.zeros((n, n))
    for j in range(Nz):
        d = [k * j + a for a in range(k + 1)]   # consecutive nodes: vertex, (mid), vertex
        M[np.ix_(d, d)] += dz * Mu
        K[np.ix_(d, d)] += Ku / dz
    inner = np.arange(1, n - 1)
    ev = sla.eigh(K[np.ix_(inner, inner)], M[np.ix_(inner, inner)], eigvals_only=True)
    return np.sort(ev)


# ----------------------------------------------------------------------------------
# dimensionless groups of the experiments
# ----------------------------------------------------------------------------------
def groups(dt, Lx, Lz, Nx, Nz, K, Nfreq=1.0e-2, k=2):
    '''
    All groups for delta = K delta_thr.  Nx is the number of columns on the FINEST level.
    '''
    eps = Lz / Lx
    m = 1.0 + (dt * Nfreq) ** 2 / 4.0
    mu = m * eps ** 2
    CM2 = max(Lx ** 2 / (4 * np.pi ** 2), m * Lz ** 2 / np.pi ** 2)
    c = CM2 / K
    dx, dz = Lx / Nx, Lz / Nz
    alpha = c / dx ** 2
    beta = c / (m * dz ** 2)
    b1 = c * (np.pi / Lz) ** 2 / m
    return dict(eps=eps, m=m, mu=mu, four_mu=4 * mu, regime='weak' if 4 * mu <= 1 else 'strong',
                c=c, gamma=c / Lx ** 2, alpha=alpha, beta=beta, b1=b1,
                mesh_aspect2_m=m * (dz / dx) ** 2, dx=dx, dz=dz,
                delta=dt / (2 * c), dtN=dt * Nfreq)


if __name__ == '__main__':
    # self test: Galerkin property and Schur complement
    for k in (1, 2):
        f = Level(16, 1.0, k)
        co = Level(8, 2.0, k)
        P = prolongation(co, f)
        c, lam, m = 37.0, 0.9, 2.5
        Af, Ac = f.matrix(c, lam, m), co.matrix(c, lam, m)
        print(f'k = {k}: Galerkin defect |P^T A_f P - A_c| = {np.abs(P.T @ Af @ P - Ac).max():.2e}')
        nu = f.nu
        S = Af[:nu, :nu] - Af[:nu, nu:] @ np.linalg.solve(Af[nu:, nu:], Af[nu:, :nu])
        print(f'k = {k}: Schur complement defect = {np.abs(S - f.schur_u(c, lam, m)).max():.2e}')
    ev = vertical_eigenvalues(20, 1.0, 2)
    print('vertical eigenvalues (k = 2, Nz = 20, Lz = 1): first three / pi^2 =', ev[:3] / np.pi ** 2,
          ' largest * dz^2 =', ev[-1] / 400.)
    print(groups(50., 3e5, 4e3, 80, 20, 1.0))
