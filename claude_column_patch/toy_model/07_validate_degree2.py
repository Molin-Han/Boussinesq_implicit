'''
Validation of the degree 2 toy model (rtcf2_model.py).  Writes 07_validate_degree2.txt.

  T1  spectrum of B A: direct dense == block Fourier
  T2  spectrum of B A against Firedrake (RTCF2, ASMStarPC / ASMVankaPC), read from
      tables/firedrake_crosscheck_spectra.npz written by 04_firedrake_crosscheck.py
  T3  Galerkin coarse operator == rediscretisation
  T4  dim ker(div) = 2 Nx (2 Nz - 1) + 1; divergence-free functions in one patch:
      star 3 (2 Nz - 1), column (2 Nz - 1) (only the ones that are bubbles in x)
  T5  multigrid FGMRES reaches the tolerance
'''
import os
import sys
import numpy as np
import rtcf2_model as r2
import rt0_model as rt

out = open('07_validate_degree2.txt', 'w')
ok_all = True


def log(*args):
    s = ' '.join(str(a) for a in args)
    print(s)
    out.write(s + '\n')


def check(name, cond, detail=''):
    global ok_all
    ok_all = ok_all and bool(cond)
    log(f'  [{"PASS" if cond else "FAIL"}] {name} {detail}')


Lx, Lz = 3.0e5, 4.0e3
pn = rt.physical_parameters(Lx, Lz, 50.0, 1.0)
ps = rt.physical_parameters(Lx, Lz, 1e5, 1.0)
cases = [
    ('unit square cells, c = 1', 8, 6, 8.0, 6.0, 1.0, 1.0),
    ('flat cells, large c', 8, 6, 8.0, 1.8, 2.5, 1.0e4),
    ('nominal physical numbers (dt = 50 s, K = 1), 16 x 10', 16, 10, Lx, Lz, pn['m'], pn['c']),
    ('strongly stratified (dt = 1e5 s, K = 1), 16 x 10', 16, 10, Lx, Lz, ps['m'], ps['c']),
]

log('T1 direct against block Fourier')
for name, nx, nz, lx, lz, m, c in cases[:3]:
    for kind in r2.PATCHES:
        ld = r2.spectrum_direct(nx, nz, lx / nx, lz / nz, m, c, kind)
        th, lam = r2.spectrum_block_fourier(nx, nz, lx / nx, lz / nz, m, c, kind)
        lb = np.sort(lam.ravel())
        err = np.max(np.abs(ld - lb) / ld)
        s = r2.spectral_summary(nx, nz, lx / nx, lz / nz, m, c, kind)
        consistent = (abs(s['lmin'] - ld[0]) < 1e-7 * ld[0] and abs(s['lmax'] - ld[-1]) < 1e-7 * ld[-1]
                      and s['ndof'] == ld.size and s['n_small'] == int(np.sum(ld < 0.1 * ld[-1])))
        check(f'{name:55s} {kind:6s}', err < 1e-6 and consistent,
              f'n = {ld.size}, lmin {ld[0]:.6e}, lmax {ld[-1]:.6f}, max rel. dev. {err:.1e}, summary consistent {consistent}')

log('\nT2 against Firedrake (degree 2)')
fn = 'tables/firedrake_crosscheck_spectra.npz'
if os.path.exists(fn):
    ref = np.load(fn)
    for name, nx, nz, lx, lz, m, c in cases:
        for kind in r2.PATCHES:
            lf = ref[f'deg2__{name}__{kind}']
            th, lam = r2.spectrum_block_fourier(nx, nz, lx / nx, lz / nz, m, c, kind)
            lb = np.sort(lam.ravel())
            same_size = lf.size == lb.size
            err = np.max(np.abs(lf - lb) / lb) if same_size else np.inf
            check(f'{name:55s} {kind:6s}', same_size and err < 1e-6,
                  f'n = {lb.size} (Firedrake {lf.size}), lmin {lb[0]:.6e} (Firedrake {lf[0]:.6e}), '
                  f'lmax {lb[-1]:.6f} (Firedrake {lf[-1]:.6f}), max rel. dev. {err:.1e}')
else:
    log('  reference file not found, run 04_firedrake_crosscheck.py first')
    ok_all = False

log('\nT3 Galerkin coarse operator equals rediscretisation')
for name, nx, nz, lx, lz, m, c in cases:
    Ac = r2.assemble(nx, nz, lx / nx, lz / nz, m, c)
    Af = r2.assemble(2 * nx, nz, lx / (2 * nx), lz / nz, m, c)
    P = r2.prolongation(nx, nz)
    G = (P.T @ Af @ P).toarray()
    err = np.abs(G - Ac.toarray()).max() / np.abs(Ac.toarray()).max()
    check(name, err < 1e-12, f'max rel. difference {err:.2e}')

log('\nT4 kernel of the divergence')
for (nx, nz, dx, dz) in [(6, 4, 1.0, 1.0), (6, 4, 3750.0, 200.0)]:
    A, M, D = r2.assemble(nx, nz, dx, dz, 1.0, 1.0, return_parts=True)
    Dd = D.toarray()
    Dd = Dd / np.abs(Dd).max(axis=0, keepdims=True)
    s = np.linalg.svd(Dd, compute_uv=False)
    dimker = Dd.shape[1] - int(np.sum(s > 1e-10 * s[0]))
    check(f'Nx={nx} Nz={nz} dx/dz={dx / dz:g}: dim ker(div)', dimker == 2 * nx * (2 * nz - 1) + 1,
          f'= {dimker}, expected 2 Nx (2 Nz - 1) + 1 = {2 * nx * (2 * nz - 1) + 1}; ndof = {Dd.shape[1]}')
    for kind, expected in [('star', 3 * (2 * nz - 1)), ('column', 2 * nz - 1)]:
        idx = r2.patch_indices(nx, nz, kind)
        Dloc = Dd[:, idx[2]]
        s = np.linalg.svd(Dloc, compute_uv=False)
        dk = Dloc.shape[1] - int(np.sum(s > 1e-10 * s[0]))
        check(f'   {kind:6s} patch size {idx.shape[1]}: dim(ker div in patch)', dk == expected, f'= {dk}, expected {expected}')
    for kind in r2.PATCHES:
        idx = r2.patch_indices(nx, nz, kind)
        basis = []
        for i in range(nx):
            Dloc = Dd[:, idx[i]]
            u, s, vt = np.linalg.svd(Dloc)
            for v in vt[np.sum(s > 1e-10 * s[0]):]:
                g = np.zeros(Dd.shape[1])
                g[idx[i]] = v
                basis.append(g)
        rank = np.linalg.matrix_rank(np.array(basis), tol=1e-10)
        expected = 2 * nx * (2 * nz - 1) if kind == 'star' else nx * (2 * nz - 1)
        check(f'   {kind:6s} sum of the local kernels', rank == expected,
              f'has dimension {rank} of dim ker(div) = {dimker} (expected {expected})')

log('\nT5 multigrid FGMRES')
name, nx, nz, lx, lz, m, c = cases[2]
for kind in r2.PATCHES:
    for cycle in ['v', 'full']:
        its, conv, res = r2.mg_fgmres_iterations(nx, nz, lx, lz, m, c, kind, cycle=cycle)
        check(f'{kind:6s} {cycle:4s}', conv and res < 1.01e-8, f'iterations {its}, true relative residual {res:.2e}')

log('\nALL TESTS PASSED' if ok_all else '\nSOME TESTS FAILED')
out.close()
sys.exit(0 if ok_all else 1)
