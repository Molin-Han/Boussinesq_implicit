'''
Validation of the toy model (rt0_model.py).  Writes 01_validate.txt.

  T1  assembly: Kronecker assembly == independent element-by-element quadrature assembly
  T2  dim ker(div) = Nx (Nz - 1) + 1 for every aspect ratio (it is a topological number)
  T3  divergence-free functions inside one patch: star Nz - 1, column 0;
      the star patches together span the kernel up to the uniform horizontal flow
  T4  Galerkin coarse operator P^T A P == rediscretised coarse operator
  T5  spectrum of B A: direct dense == block Fourier in x == 2 x 2 symbols
  T6  ASM.apply == dense B
  T7  smoother GMRES: result == dense least-squares minimiser over the Krylov space
  T8  FGMRES: true residual below tolerance, one iteration with the exact inverse
  T9  closed forms: column patch has the Nz-fold eigenvalue 2 / (5/6 + c/dx^2) at theta = 0,
      bounds lambda_max <= 2 (star) and <= 3 (column)
'''
import sys
import numpy as np
import scipy.sparse as sp
import scipy.linalg as sla
import rt0_model as rt

out = open('01_validate.txt', 'w')


def log(*args):
    s = ' '.join(str(a) for a in args)
    print(s)
    out.write(s + '\n')


ok_all = True


def check(name, cond, detail=''):
    global ok_all
    ok_all = ok_all and bool(cond)
    log(f'  [{"PASS" if cond else "FAIL"}] {name} {detail}')


# nominal physical numbers
Lx, Lz = 3.0e5, 4.0e3
pp = rt.physical_parameters(Lx, Lz, dt=50.0, K=1.0)
log('nominal parameters (dt = 50 s, K = 1):', pp)

cases = [
    dict(name='unit', Nx=8, Nz=6, dx=1.0, dz=1.0, m=1.0, c=1.0),
    dict(name='unit, large c', Nx=8, Nz=6, dx=1.0, dz=0.3, m=2.5, c=1.0e4),
    dict(name='nominal physical, coarse mesh', Nx=8, Nz=6, dx=Lx / 8, dz=Lz / 6, m=pp['m'], c=pp['c']),
    dict(name='nominal physical, 16 x 10', Nx=16, Nz=10, dx=Lx / 16, dz=Lz / 10, m=pp['m'], c=pp['c']),
    dict(name='strongly stratified (dt = 1e5), 16 x 10', Nx=16, Nz=10, dx=Lx / 16, dz=Lz / 10,
         m=rt.physical_parameters(Lx, Lz, 1e5, 1.0)['m'], c=rt.physical_parameters(Lx, Lz, 1e5, 1.0)['c']),
]

log('\nT1 assembly against quadrature')
for cs in cases[:3]:
    A = rt.assemble(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c']).toarray()
    Aq = rt.assemble_by_quadrature(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'])
    err = np.abs(A - Aq).max() / np.abs(A).max()
    check(cs['name'], err < 1e-13, f'max rel. difference {err:.2e}')
    check(cs['name'] + ' symmetric positive definite', np.allclose(A, A.T) and np.linalg.eigvalsh(A).min() > 0)

log('\nT2 dimension of ker(div) for several aspect ratios')
for (Nx, Nz, dx, dz) in [(8, 6, 1.0, 1.0), (8, 6, 1.0, 1e-3), (8, 6, 1e-3, 1.0), (12, 5, 3750.0, 200.0)]:
    A, M, D = rt.assemble(Nx, Nz, dx, dz, 1.0, 1.0, return_parts=True)
    Dd = D.toarray()
    # scale the columns so that the rank decision does not depend on the units
    s = np.linalg.svd(Dd / np.abs(Dd).max(axis=0, keepdims=True), compute_uv=False)
    rank = int(np.sum(s > 1e-10 * s[0]))
    dimker = Dd.shape[1] - rank
    check(f'Nx={Nx} Nz={Nz} dx/dz={dx / dz:g}', dimker == Nx * (Nz - 1) + 1,
          f'dim ker = {dimker}, expected {Nx * (Nz - 1) + 1}; rank(div) = {rank} = Nx*Nz - 1 = {Nx * Nz - 1}')

log('\nT3 divergence-free functions supported in one patch')
for (Nx, Nz, dx, dz) in [(8, 6, 1.0, 1.0), (12, 5, 3750.0, 200.0)]:
    A, M, D = rt.assemble(Nx, Nz, dx, dz, 1.0, 1.0, return_parts=True)
    Dd = D.toarray()
    Dd = Dd / np.abs(Dd).max(axis=0, keepdims=True)
    for kind, expected in [('star', Nz - 1), ('column', 0)]:
        idx = rt.patch_indices(Nx, Nz, kind)
        Dloc = Dd[:, idx[2]]
        s = np.linalg.svd(Dloc, compute_uv=False)
        dimker = Dloc.shape[1] - int(np.sum(s > 1e-10 * s[0]))
        check(f'{kind:6s} Nx={Nx} Nz={Nz} dx/dz={dx / dz:g}: patch size {idx.shape[1]}', dimker == expected,
              f'dim(ker div in patch) = {dimker}, expected {expected}')
    # span of the local kernels of all star patches
    idx = rt.patch_indices(Nx, Nz, 'star')
    basis = []
    for i in range(Nx):
        Dloc = Dd[:, idx[i]]
        u, s, vt = np.linalg.svd(Dloc)
        nullvecs = vt[np.sum(s > 1e-10 * s[0]):]
        for v in nullvecs:
            g = np.zeros(Dd.shape[1])
            g[idx[i]] = v
            basis.append(g)
    basis = np.array(basis)
    rank = np.linalg.matrix_rank(basis, tol=1e-10)
    check(f'star   Nx={Nx} Nz={Nz}: sum of the local kernels', rank == Nx * (Nz - 1),
          f'has dimension {rank} = dim ker(div) - 1 = {Nx * (Nz - 1)} (the missing one is the uniform flow (1, 0))')

log('\nT4 Galerkin coarse operator equals rediscretisation')
for cs in cases:
    Nxc, Nz = cs['Nx'], cs['Nz']
    Ac = rt.assemble(Nxc, Nz, cs['dx'], cs['dz'], cs['m'], cs['c'])
    Af = rt.assemble(2 * Nxc, Nz, cs['dx'] / 2, cs['dz'], cs['m'], cs['c'])
    P = rt.prolongation(Nxc, Nz)
    G = (P.T @ Af @ P).toarray()
    err = np.abs(G - Ac.toarray()).max() / np.abs(Ac.toarray()).max()
    check(cs['name'], err < 1e-12, f'max rel. difference {err:.2e}')

log('\nT5 spectrum of B A: direct, block Fourier, 2 x 2 symbols')
for cs in cases:
    for kind in rt.PATCHES:
        args = (cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'], kind)
        ld = rt.spectrum_direct(*args)
        lb = rt.spectrum_block_fourier(*args)
        l2 = rt.spectrum_2x2(*args)
        e1 = np.max(np.abs(ld - lb) / ld)
        e2 = np.max(np.abs(ld - l2) / ld)
        check(f'{cs["name"]:42s} {kind:6s}', e1 < 1e-7 and e2 < 1e-7,
              f'n = {ld.size}, lmin = {ld[0]:.6e}, lmax = {ld[-1]:.6f}, cond = {ld[-1] / ld[0]:.4e}; '
              f'max rel. dev. block Fourier {e1:.1e}, 2x2 symbols {e2:.1e}')

log('\nT6 ASM.apply against the dense matrix')
rng = np.random.default_rng(1)
for cs in cases[2:4]:
    A = rt.assemble(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'])
    for kind in rt.PATCHES:
        asm = rt.ASM(A, cs['Nx'], cs['Nz'], kind, check_translation_invariance=True)
        r = rng.standard_normal(A.shape[0])
        y1 = asm.apply(r)
        y2 = asm.dense() @ r
        err = np.linalg.norm(y1 - y2) / np.linalg.norm(y2)
        check(f'{cs["name"]} {kind}', err < 1e-12, f'rel. difference {err:.1e}')

log('\nT7 smoother: GMRES(6) against the dense least-squares minimiser')
for cs in cases[2:4]:
    A = rt.assemble(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'])
    for kind in rt.PATCHES:
        asm = rt.ASM(A, cs['Nx'], cs['Nz'], kind)
        Bd = asm.dense()
        Ad = A.toarray()
        b = rng.standard_normal(A.shape[0])
        x0 = rng.standard_normal(A.shape[0]) * 1e-7
        x = rt.gmres_fixed(lambda v: A @ v, asm.apply, b, x0, 6)
        r0 = Bd @ (b - Ad @ x0)
        Kr = [r0]
        for _ in range(5):
            Kr.append(Bd @ (Ad @ Kr[-1]))
        Q, _ = np.linalg.qr(np.array(Kr).T)
        y, *_ = np.linalg.lstsq(Bd @ Ad @ Q, r0, rcond=None)
        xref = x0 + Q @ y
        n1 = np.linalg.norm(Bd @ (b - Ad @ x))
        n2 = np.linalg.norm(Bd @ (b - Ad @ xref))
        check(f'{cs["name"]} {kind}', abs(n1 - n2) / n2 < 1e-6,
              f'preconditioned residual norms {n1:.8e} (gmres_fixed) {n2:.8e} (least squares)')

log('\nT8 FGMRES')
cs = cases[3]
A = rt.assemble(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'])
Ad = A.toarray()
b = rng.standard_normal(A.shape[0])
lu = sla.lu_factor(Ad)
x, its, hist, conv = rt.fgmres(lambda v: A @ v, lambda v: sla.lu_solve(lu, v), b, rtol=1e-8)
check('exact inverse as preconditioner', conv and its == 1, f'iterations {its}')
for kind in rt.PATCHES:
    asm = rt.ASM(A, cs['Nx'], cs['Nz'], kind)
    x, its, hist, conv = rt.fgmres(lambda v: A @ v, asm.apply, b, rtol=1e-8, maxit=2000, restart=400)
    res = np.linalg.norm(b - A @ x) / np.linalg.norm(b)
    check(f'one-level ASM {kind}', conv and res < 1.01e-8, f'iterations {its}, true relative residual {res:.2e}')
for kind in rt.PATCHES:
    for cycle in ['v', 'full']:
        its, conv, res = rt.mg_fgmres_iterations(16, 10, Lx, Lz, cs['m'], cs['c'], kind, nlevels=3, cycle=cycle)
        check(f'multigrid {cycle:4s} {kind:6s}', conv and res < 1.01e-8, f'iterations {its}, true relative residual {res:.2e}')

log('\nT9 closed forms')
for cs in cases:
    gamma = cs['c'] / cs['dx'] ** 2
    th, lam1, lam0 = rt.symbols_2x2(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'], 'column')
    target = 2.0 / (5.0 / 6.0 + gamma)
    lam_theta0 = np.sort(np.concatenate([lam1[0].ravel(), lam0[:1]]))
    n_match = int(np.sum(np.abs(lam_theta0 - target) < 1e-8 * target))
    check(f'{cs["name"]:42s} column: eigenvalue 2/(5/6 + c/dx^2) = {target:.6e} at theta = 0',
          n_match == cs['Nz'], f'multiplicity found {n_match}, expected Nz = {cs["Nz"]}')
    for kind, bound in [('star', 2.0), ('column', 3.0)]:
        l = rt.spectrum_2x2(cs['Nx'], cs['Nz'], cs['dx'], cs['dz'], cs['m'], cs['c'], kind)
        check(f'{cs["name"]:42s} {kind:6s}: lambda_max <= {bound:g}', l[-1] <= bound * (1 + 1e-9), f'lambda_max = {l[-1]:.9f}')

log('\nALL TESTS PASSED' if ok_all else '\nSOME TESTS FAILED')
out.close()
sys.exit(0 if ok_all else 1)
