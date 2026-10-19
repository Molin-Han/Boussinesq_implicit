'''
Reduced (modal) model of the linear solver of Nonlinear_Boussinesq_Irk_SC.py.

What is modelled
----------------
The Jacobian of the implicit-midpoint stage problem linearised about the state of rest
with f = 0, in the basis of vertical modes (see modal1d.py; the reduction is exact and is
verified against Firedrake in verify_modal_reduction_firedrake.py).  For a vertical mode
with eigenvalue lambda the unknowns are functions of x only,

    U in CG_2(x)   (horizontal velocity, vertical structure psi_n in DG_1(z)),
    W in DG_1(x)   (vertical velocity,   vertical structure chi_n in CG_2(z)),
    B in DG_1(x)   (buoyancy,            vertical structure chi_n),
    P in DG_1(x)   (pressure,            vertical structure psi_n),

and the stage-derivative Jacobian is (h = dt/2)

    [ Mu            0              0        -h C^T       ] [U]
    [ 0             Mw            -h Mw     -h lam Mw    ] [W]
    [ 0             h N^2 Mw       Mw        0           ] [B]
    [ h C           h lam Mw       0         0           ] [P]

The preconditioner is the solver's: auxiliary operator with p -> p - div(u)/delta and
delta p phi dx (class HDivSchurPC), i.e. block upper triangular

    [ A11   -h G ]        A11 = field-1 operator with c (div, div), c = dt / (2 delta),
    [ 0   delta h Mp ]

applied as: pressure mass solve, then ONE application of full multigrid (PETSc PCMG type
FULL, V cycle, x-coarsening, rediscretised coarse operators, coarse LU) on A11 with
smoother = GMRES(6), left preconditioned by the additive Schwarz method with star or vanka
patches (exact patch solves).  The outer iteration is right-preconditioned FGMRES(30),
rtol 1e-6, atol 1e-9, maxit 100, as in the script.

What is NOT modelled
--------------------
* the nonlinear (advective) terms of the Jacobian and of the residual,
* the out-of-plane velocity u_y (decoupled for f = 0),
* the Euclidean norms of the production vectors are reproduced only up to a few per cent:
  the basis of the emulator is scaled like the basis of the solver (contravariant Piola map:
  the velocity dofs are fluxes, |phi_u| ~ 1/dz, |phi_w| ~ 1/dx; DG_1 Lagrange basis at the
  Gauss points, mass (h/2) I), i.e. the modal basis functions are multiplied by
      s_U = 1/sqrt(2 dz),  s_W = sqrt(dz/2)/dx,  s_B = s_P = sqrt(dz/2);
  this is exact for the DG_1(z) fields (u, p) and approximate (5 %) for the CG_2(z) fields.
* parallel decomposition (2 ranks) -- the solver's ASM is the same in parallel.

Data: second time step of the solver's test (rest state + buoyancy perturbation
3e-4 sin(pi z / L_z) / (1 + (x - L_x/2)^2 / a^2), a = 5 km), first step solved exactly.
'''
import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import modal1d

SQ2 = np.sqrt(2.0)


# ==================================================================================
class ModeLevel:
    '''Matrices of one vertical mode on one level (N periodic cells of width dx, degree 2).'''

    def __init__(self, N, dx, dz=2.0, scaled=True):
        self.N, self.dx, self.dz = N, dx, dz
        if scaled:
            self.sU, self.sW = 1.0 / np.sqrt(2.0 * dz), np.sqrt(dz / 2.0) / dx
            self.sB = self.sP = np.sqrt(dz / 2.0)
        else:
            self.sU = self.sW = self.sB = self.sP = 1.0
        L = modal1d.Level(N, dx, 2)
        self.L = L
        self.nu = 2 * N            # U: vertices 0..N-1, midpoints N..2N-1
        self.nd = 2 * N            # each DG_1 field
        self.Mu = sp.csr_matrix(L.Mu)
        self.Ku = sp.csr_matrix(L.Ku)
        self.C = sp.csr_matrix(L.C / SQ2)             # DG basis with mass (dx/2) I
        self.Md = sp.identity(self.nd, format='csr') * (dx / 2.0)
        self.n1 = self.nu + 2 * self.nd               # field 1 = (U, W, B)
        self.nall = self.n1 + self.nd                 # + P
        self.D1 = sp.diags(np.concatenate([np.full(self.nu, self.sU), np.full(self.nd, self.sW),
                                           np.full(self.nd, self.sB)]), format='csr')
        self.Dall = sp.diags(np.concatenate([np.full(self.nu, self.sU), np.full(self.nd, self.sW),
                                             np.full(self.nd, self.sB), np.full(self.nd, self.sP)]),
                             format='csr')

    # index helpers (field-1 numbering)
    def iU(self):
        return np.arange(self.nu)

    def iW(self):
        return self.nu + np.arange(self.nd)

    def iB(self):
        return self.nu + self.nd + np.arange(self.nd)

    def cell_dg(self, i):
        return [2 * i, 2 * i + 1]

    def A11(self, c, lam, dt, N2):
        h = dt / 2.0
        Z = None
        A = sp.bmat([[self.Mu + c * self.Ku, c * lam * self.C.T, Z],
                     [c * lam * self.C, (1.0 + c * lam ** 2) * self.Md, -h * self.Md],
                     [Z, h * N2 * self.Md, self.Md]], format='csr')
        return (self.D1 @ A @ self.D1).tocsr()

    def A11_plain(self, dt, N2):
        '''field-1 block of the true Jacobian (no grad-div).'''
        h = dt / 2.0
        Z = sp.csr_matrix((self.nd, self.nu))
        A = sp.bmat([[self.Mu, Z.T, Z.T],
                     [Z, self.Md, -h * self.Md],
                     [Z, h * N2 * self.Md, self.Md]], format='csr')
        return (self.D1 @ A @ self.D1).tocsr()

    def G(self, lam):
        '''G P: the pressure-gradient coupling; rows (U, W, B), column P:  [C^T; lam Md; 0].'''
        G = sp.vstack([self.C.T, lam * self.Md, sp.csr_matrix((self.nd, self.nd))], format='csr')
        return (self.D1 @ G * self.sP).tocsr()

    def mass_all(self):
        M = sp.block_diag([self.Mu, self.Md, self.Md, sp.csr_matrix((self.nd, self.nd))], format='csr')
        return (self.Dall @ M @ self.Dall).tocsr()

    def patches(self, kind):
        N, nu, nd = self.N, self.nu, self.nd
        out = []
        for i in range(N):
            if kind == 'star':
                d = [i]
                for cell in ((i - 1) % N, i):
                    d.append(N + cell)
                    d += [nu + j for j in self.cell_dg(cell)]
                    d += [nu + nd + j for j in self.cell_dg(cell)]
            elif kind == 'vanka':
                d = [i, (i + 1) % N, N + i]
                d += [nu + j for j in self.cell_dg(i)]
                d += [nu + nd + j for j in self.cell_dg(i)]
            else:
                raise ValueError(kind)
            out.append(np.array(sorted(set(d))))
        return out


def prolong_field1(coarse, fine):
    '''Natural embedding for (U, W, B).'''
    Pfull = modal1d.prolongation(coarse.L, fine.L, with_w=True)   # (U, W) with orthonormal DG basis
    nuf, nuc = fine.nu, coarse.nu
    Pu = Pfull[:nuf, :nuc]
    # DG: both bases are scaled by the same constant factor relative to the orthonormal one
    # per level: psi' = psi / sqrt(2) on both levels, but 'orthonormal' means mass dx I, and
    # dx_c = 2 dx_f, so the coarse function sum_c Wc psi_c' equals sum_f (P Wc) psi_f' with the
    # same P as for the basis with mass dx I.
    Pd = Pfull[nuf:, nuc:]
    P = sp.block_diag([sp.csr_matrix(Pu), sp.csr_matrix(Pd), sp.csr_matrix(Pd)], format='csr')
    Dfi = sp.diags(1.0 / fine.D1.diagonal())
    return (Dfi @ P @ coarse.D1).tocsr()


# ==================================================================================
class ASM:
    def __init__(self, A, patches):
        self.patches = patches
        Ad = A.toarray() if sp.issparse(A) else A
        self.inv = [np.linalg.inv(Ad[np.ix_(p, p)]) for p in patches]

    def apply(self, r):
        y = np.zeros_like(r)
        for p, Ai in zip(self.patches, self.inv):
            y[p] += Ai @ r[p]
        return y


def gmres_left(A, B, b, x0, maxit=6, rtol=0.0, atol=0.0):
    '''
    Left-preconditioned GMRES as PETSc's KSPGMRES: minimises ||B (b - A x)||_2 over the Krylov
    space of B A; the convergence test uses the preconditioned residual norm relative to ||B b||.
    A, B: callables.  Returns x.
    '''
    r = B(b - A(x0)) if np.any(x0) else B(b)
    beta = np.linalg.norm(r)
    ref = np.linalg.norm(B(b)) if np.any(x0) else beta
    if beta <= max(rtol * ref, atol) or beta == 0.0:
        return x0.copy()
    V = [r / beta]
    H = np.zeros((maxit + 1, maxit))
    y = None
    for j in range(maxit):
        w = B(A(V[j]))
        for i in range(j + 1):
            H[i, j] = np.dot(V[i], w)
            w = w - H[i, j] * V[i]
        H[j + 1, j] = np.linalg.norm(w)
        e1 = np.zeros(j + 2)
        e1[0] = beta
        y, res, _, _ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
        rn = np.linalg.norm(H[:j + 2, :j + 1] @ y - e1)
        if rn <= max(rtol * ref, atol) or H[j + 1, j] == 0.0:
            break
        V.append(w / H[j + 1, j])
    x = x0.copy()
    for i in range(len(y)):
        x = x + y[i] * V[i]
    return x


def fgmres_right(A, Minv, b, rtol=1e-6, atol=1e-9, maxit=100, restart=30, norm_scale=1.0):
    '''
    Right-preconditioned flexible GMRES with zero initial guess.
    The convergence test is on norm_scale * ||b - A x||_2 (unpreconditioned), as in the solver.
    Returns x, list of scaled residual norms (iteration 0 .. n), reason.
    '''
    x = np.zeros_like(b)
    r = b.copy()
    rn0 = np.linalg.norm(r) * norm_scale
    hist = [rn0]
    target = max(rtol * rn0, atol)
    if rn0 <= atol:
        return x, hist, 'atol'
    it = 0
    while it < maxit:
        beta = np.linalg.norm(r)
        V = [r / beta]
        Z = []
        H = np.zeros((restart + 1, restart))
        y = None
        for j in range(restart):
            z = Minv(V[j])
            Z.append(z)
            w = A(z)
            for i in range(j + 1):
                H[i, j] = np.dot(V[i], w)
                w = w - H[i, j] * V[i]
            H[j + 1, j] = np.linalg.norm(w)
            e1 = np.zeros(j + 2)
            e1[0] = beta
            y, _, _, _ = np.linalg.lstsq(H[:j + 2, :j + 1], e1, rcond=None)
            rn = np.linalg.norm(H[:j + 2, :j + 1] @ y - e1) * norm_scale
            it += 1
            hist.append(rn)
            if rn <= target or it >= maxit:
                break
            V.append(w / H[j + 1, j])
        for i in range(len(y)):
            x = x + y[i] * Z[i]
        r = b - A(x)
        if hist[-1] <= target:
            return x, hist, ('rtol' if rtol * rn0 >= atol else 'atol')
    return x, hist, 'maxit'


# ==================================================================================
class ModalSolver:
    '''
    The solver for a set of vertical modes.  Vectors are the concatenation over the modes.
    '''

    def __init__(self, dt, Lx, Lz, nx_coarse, Nz, K, patch, refinement=2, Nfreq=1.0e-2,
                 modes=(1,), smoother_its=6, delta=None, scaled=True):
        self.dt, self.Lx, self.Lz, self.Nz, self.K = dt, Lx, Lz, Nz, K
        self.N2 = Nfreq ** 2
        self.patch = patch
        self.nlev = refinement + 1
        Nx = nx_coarse * 2 ** refinement
        self.Nx = Nx
        g = modal1d.groups(dt, Lx, Lz, Nx, Nz, K, Nfreq)
        self.groups = g
        self.delta = g['delta'] if delta is None else delta
        self.c = dt / (2.0 * self.delta)
        self.m = g['m']
        lam2 = modal1d.vertical_eigenvalues(Nz, Lz, 2)
        self.modes = list(modes)
        self.lams = [np.sqrt(lam2[n - 1]) for n in self.modes]
        self.smoother_its = smoother_its
        self.levels = [ModeLevel(nx_coarse * 2 ** l, Lx / (nx_coarse * 2 ** l), dz=Lz / Nz, scaled=scaled)
                       for l in range(self.nlev)]
        self.P = [None] + [prolong_field1(self.levels[l - 1], self.levels[l]) for l in range(1, self.nlev)]
        fine = self.levels[-1]
        self.fine = fine
        h = dt / 2.0
        # per mode operators
        self.A11 = []      # [mode][level]
        self.asm = []
        self.lu0 = []
        for lam in self.lams:
            Al = [lev.A11(self.c, lam, dt, self.N2) for lev in self.levels]
            self.A11.append(Al)
            self.asm.append([None] + [ASM(Al[l], self.levels[l].patches(patch)) for l in range(1, self.nlev)])
            self.lu0.append(spla.splu(sp.csc_matrix(Al[0])))
        # full Jacobian and the pieces of the preconditioner on the finest level
        self.J = []
        self.Gf = []
        for lam in self.lams:
            G = fine.G(lam)
            A0 = fine.A11_plain(dt, self.N2)
            J = sp.bmat([[A0, -h * G], [h * G.T, None]], format='csr')
            # G.T has rows P and columns (U, W, B): (C U + lam Md W)
            self.J.append(J)
            self.Gf.append(G)
        self.n1 = fine.n1
        self.nall = fine.nall
        self.nmodes = len(self.modes)
        self.vertical_norm_factor = 1.0 if scaled else np.sqrt((Lz / Nz) / 2.0)

    # ---------------------------------------------------------------- vectors
    def split(self, v, size):
        return [v[i * size:(i + 1) * size] for i in range(self.nmodes)]

    # ---------------------------------------------------------------- operators
    def J_apply(self, v):
        return np.concatenate([self.J[i] @ vi for i, vi in enumerate(self.split(v, self.nall))])

    def A11_apply(self, l, v):
        n = self.levels[l].n1
        return np.concatenate([self.A11[i][l] @ v[i * n:(i + 1) * n] for i in range(self.nmodes)])

    def asm_apply(self, l, v):
        n = self.levels[l].n1
        return np.concatenate([self.asm[i][l].apply(v[i * n:(i + 1) * n]) for i in range(self.nmodes)])

    def coarse_solve(self, v):
        n = self.levels[0].n1
        return np.concatenate([self.lu0[i].solve(v[i * n:(i + 1) * n]) for i in range(self.nmodes)])

    def prolong(self, l, v):
        '''from level l-1 to level l'''
        n = self.levels[l - 1].n1
        return np.concatenate([self.P[l] @ v[i * n:(i + 1) * n] for i in range(self.nmodes)])

    def restrict(self, l, v):
        '''from level l to level l-1'''
        n = self.levels[l].n1
        return np.concatenate([self.P[l].T @ v[i * n:(i + 1) * n] for i in range(self.nmodes)])

    # ---------------------------------------------------------------- multigrid
    def smooth(self, l, b, x):
        return gmres_left(lambda v: self.A11_apply(l, v), lambda v: self.asm_apply(l, v), b, x,
                          maxit=self.smoother_its)

    def vcycle(self, l, b, x):
        if l == 0:
            return self.coarse_solve(b)
        x = self.smooth(l, b, x)
        r = b - self.A11_apply(l, x)
        bc = self.restrict(l, r)
        xc = self.vcycle(l - 1, bc, np.zeros_like(bc))
        x = x + self.prolong(l, xc)
        x = self.smooth(l, b, x)
        return x

    def full_mg(self, b):
        '''PETSc PCMG type FULL with V cycles.'''
        bs = [None] * self.nlev
        bs[-1] = b
        for l in range(self.nlev - 1, 0, -1):
            bs[l - 1] = self.restrict(l, bs[l])
        x = self.vcycle(0, bs[0], None)
        for l in range(1, self.nlev):
            x = self.prolong(l, x)
            x = self.vcycle(l, bs[l], x)
        return x

    def exact_field1(self, b):
        if not hasattr(self, '_luf'):
            self._luf = [spla.splu(sp.csc_matrix(self.A11[i][-1])) for i in range(self.nmodes)]
        n = self.n1
        return np.concatenate([self._luf[i].solve(b[i * n:(i + 1) * n]) for i in range(self.nmodes)])

    # ---------------------------------------------------------------- preconditioner
    def precond(self, r, exact=False):
        h = self.dt / 2.0
        out = []
        r1s, xps = [], []
        for i, ri in enumerate(self.split(r, self.nall)):
            r1, rp = ri[:self.n1], ri[self.n1:]
            xp = rp / (self.delta * h * self.fine.dx / 2.0 * self.fine.sP ** 2)
            r1s.append(r1 + h * (self.Gf[i] @ xp))
            xps.append(xp)
        rhs = np.concatenate(r1s)
        x1 = self.exact_field1(rhs) if exact else self.full_mg(rhs)
        for i in range(self.nmodes):
            out.append(np.concatenate([x1[i * self.n1:(i + 1) * self.n1], xps[i]]))
        return np.concatenate(out)

    # ---------------------------------------------------------------- data
    def initial_buoyancy(self, amplitude=3.0e-4, a=5000.0):
        '''
        L2 projection of amplitude / (1 + (x - Lx/2)^2 / a^2) onto DG_1(x) in the emulator's basis
        (psi_0 = 1/sqrt2, psi_1 = sqrt(3/2)(2 xi - 1) on each cell, mass (dx/2) I),
        times the vertical coefficients beta_n = (sin(pi z / Lz), chi_n).
        '''
        fine = self.fine
        xg, wg = np.polynomial.legendre.leggauss(12)
        xi = 0.5 * (xg + 1.0)
        wq = 0.5 * wg
        psi = np.array([np.ones_like(xi) / SQ2, np.sqrt(1.5) * (2 * xi - 1)])
        Bx = np.zeros(fine.nd)
        for i in range(fine.N):
            x = (i + xi) * fine.dx
            f = amplitude / (1.0 + (x - self.Lx / 2.0) ** 2 / a ** 2)
            for q in range(2):
                Bx[2 * i + q] = fine.dx * np.sum(wq * psi[q] * f) / (fine.dx / 2.0)
        betas = vertical_projection_of_sine(self.Nz, self.Lz)
        return [Bx * betas[n - 1] for n in self.modes]

    def second_step_problem(self):
        '''
        Returns (rhs, info): rhs = -F at the initial guess of the first Newton step of the SECOND
        time step.  State vector per mode (U, W, B, P); the spatial operator is L = (J - Mass)/h.
        '''
        h = self.dt / 2.0
        B0 = self.initial_buoyancy()
        rhs_all = []
        info = dict(step1_residual=0.0)
        r1sq = 0.0
        for i in range(self.nmodes):
            J = self.J[i]
            fine = self.fine
            mass = fine.mass_all()
            Lop = (J - mass) / h
            U0 = np.zeros(self.nall)
            U0[fine.nu + fine.nd:fine.nu + 2 * fine.nd] = B0[i] / fine.sB     # coefficients in the scaled basis
            F0 = Lop @ U0                      # residual of step 1 at k = 0
            r1sq += np.dot(F0, F0)
            # pressure nullspace does not occur for n >= 1
            k1 = spla.spsolve(sp.csc_matrix(J), -F0)
            U1 = U0 + self.dt * k1
            F1 = J @ k1 + Lop @ U1             # residual of step 2 at the initial guess k = k1
            rhs_all.append(-F1)
        info['step1_residual'] = np.sqrt(r1sq) * self.vertical_norm_factor
        rhs = np.concatenate(rhs_all)
        info['step2_residual'] = np.linalg.norm(rhs) * self.vertical_norm_factor
        return rhs, info

    # ---------------------------------------------------------------- the two linear solves
    def run(self, exact_inner=False, rtol=1e-6, atol=1e-9, maxit=100, snes_atol=1e-9, snes_rtol=1e-8,
            max_newton=10, verbose=False):
        rhs, info = self.second_step_problem()
        A = self.J_apply

        def Minv(v):
            return self.precond(v, exact=exact_inner)
        counts, hists = [], []
        r = rhs.copy()
        f0 = np.linalg.norm(r) * self.vertical_norm_factor
        for newton in range(max_newton):
            fn = np.linalg.norm(r) * self.vertical_norm_factor
            if newton > 0 and (fn <= snes_atol or fn <= snes_rtol * f0):
                break
            x, hist, reason = fgmres_right(A, Minv, r, rtol=rtol, atol=atol, maxit=maxit,
                                           norm_scale=self.vertical_norm_factor)
            counts.append(len(hist) - 1)
            hists.append(hist)
            r = r - A(x)
            if verbose:
                print(f'   Newton step {newton + 1}: {len(hist) - 1} its ({reason}), '
                      f'residual {hist[0]:.3e} -> {hist[-1]:.3e}')
            if reason == 'maxit':
                break
        info.update(counts=counts, cumulative=list(np.cumsum(counts)), histories=hists)
        return info


def vertical_projection_of_sine(Nz, Lz):
    '''beta_n = int chi_n(z) sin(pi z / Lz) dz for the M_z-orthonormal discrete modes chi_n of CG_2.'''
    dz = Lz / Nz
    Mu, Ku, _, _ = modal1d.reference(2)
    n = 2 * Nz + 1
    M = np.zeros((n, n))
    K = np.zeros((n, n))
    xg, wg = np.polynomial.legendre.leggauss(10)
    xi = 0.5 * (xg + 1)
    wq = 0.5 * wg
    phi = np.array([(1 - xi) * (1 - 2 * xi), 4 * xi * (1 - xi), xi * (2 * xi - 1)])
    F = np.zeros(n)
    for j in range(Nz):
        d = [2 * j, 2 * j + 1, 2 * j + 2]
        M[np.ix_(d, d)] += dz * Mu
        K[np.ix_(d, d)] += Ku / dz
        s = np.sin(np.pi * (j + xi) * dz / Lz)
        for a in range(3):
            F[2 * j + a] += dz * np.sum(wq * phi[a] * s)
    inner = np.arange(1, n - 1)
    ev, X = sla.eigh(K[np.ix_(inner, inner)], M[np.ix_(inner, inner)])
    beta = X.T @ F[inner]
    beta[0] = abs(beta[0])
    return beta


if __name__ == '__main__':
    import time
    for patch in ('star', 'vanka'):
        t0 = time.time()
        s = ModalSolver(dt=50.0, Lx=3.0e5, Lz=4.0e3, nx_coarse=20, Nz=20, K=1.0, patch=patch)
        info = s.run(verbose=True)
        print(patch, 'step-1 residual', info['step1_residual'], 'step-2 residual', info['step2_residual'],
              'cumulative KSP', [0] + info['cumulative'], f'{time.time() - t0:.1f} s')
        for hst in info['histories']:
            print('   rel. residuals:', ' '.join(f'{v / hst[0]:.2e}' for v in hst))
    s = ModalSolver(dt=50.0, Lx=3.0e5, Lz=4.0e3, nx_coarse=20, Nz=20, K=1.0, patch='star')
    info = s.run(exact_inner=True, verbose=True)
    print('exact field-1 solve: cumulative KSP', [0] + info['cumulative'])
