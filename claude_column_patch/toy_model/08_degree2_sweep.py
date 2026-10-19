'''
Degree 2 toy model (rtcf2_model.py): spectrum of the one-level preconditioned operator and
FGMRES iteration counts with the multigrid preconditioner (V cycle and PETSc "full" cycle),
for both patches, on reduced sweeps around the nominal configuration

    N_x = 80 (finest level), N_z = 20, L_x = 300 km, L_z = 4 km, N = 1e-2 1/s.

usage: python 08_degree2_sweep.py <group>
    dtK   : dt in {1, ..., 1e5} x K in {0.1, 1, 10, 100, 1000}
    LzA   : L_z in {12, ..., 0.25} km x K in {1, 10, 100} at dt = 50 s
    LzB   : L_z in {12, ..., 0.25} km x K in {1, 10, 100} at dt = 5000 s
    LzC   : L_z in {12, ..., 0.25} km x K in {1, 10, 100} at dt = 1e5 s
    nx    : N_x in {20, 40, 80, 160, 320} x K in {1, 10, 100} at dt = 50 s
    nz    : N_z in {10, 20, 40} x K in {1, 10, 100} at dt = 50 s

Output: tables/degree2_<group>.csv
'''
import csv
import os
import sys
import time
import itertools
import numpy as np
import rt0_model as rt
import rtcf2_model as r2

LX = 3.0e5
NFREQ = 1.0e-2
LZS = [12e3, 8e3, 4e3, 2e3, 1e3, 0.5e3, 0.25e3]
DTS = [1., 10., 50., 100., 500., 1e3, 5e3, 1e4, 5e4, 1e5]
KS = [0.1, 1., 10., 100., 1000.]
K3 = [1., 10., 100.]

FIELDS = ['group', 'degree', 'Nx', 'Nz', 'Lx', 'Lz', 'dt', 'K', 'eps', 'm', 'c', 'regime', 'four_eps2_m',
          'dx', 'dz', 'gamma', 'eta', 'patch', 'patch_size', 'ndof',
          'lmin', 'lmax', 'cond', 'n_small', 'frac_small', 'lmin_nh', 'cond_nh',
          'lmin_hf', 'cond_hf', 'n_small_hf', 'frac_small_hf',
          'its_v', 'conv_v', 'its_full', 'conv_full', 'seconds']


def configs(group):
    if group == 'dtK':
        return [(80, 20, 4e3, dt, K) for dt, K in itertools.product(DTS, KS)]
    if group in ('LzA', 'LzB', 'LzC'):
        dt = dict(LzA=50., LzB=5000., LzC=1e5)[group]
        return [(80, 20, Lz, dt, K) for Lz, K in itertools.product(LZS, K3)]
    if group == 'nx':
        return [(Nx, 20, 4e3, 50., K) for Nx, K in itertools.product([20, 40, 80, 160, 320], K3)]
    if group == 'nz':
        return [(80, Nz, 4e3, 50., K) for Nz, K in itertools.product([10, 20, 40], K3)]
    raise ValueError(group)


def run(group, maxit=300):
    os.makedirs('tables', exist_ok=True)
    fname = f'tables/degree2_{group}.csv'
    t00 = time.time()
    with open(fname, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for (Nx, Nz, Lz, dt, K) in configs(group):
            pp = rt.physical_parameters(LX, Lz, dt, K, N=NFREQ)
            dx, dz = LX / Nx, Lz / Nz
            for kind in r2.PATCHES:
                t0 = time.time()
                s = r2.spectral_summary(Nx, Nz, dx, dz, pp['m'], pp['c'], kind)
                iv, cv, _ = r2.mg_fgmres_iterations(Nx, Nz, LX, Lz, pp['m'], pp['c'], kind, cycle='v', maxit=maxit)
                ifu, cf, _ = r2.mg_fgmres_iterations(Nx, Nz, LX, Lz, pp['m'], pp['c'], kind, cycle='full', maxit=maxit)
                row = dict(group=group, degree=2, Nx=Nx, Nz=Nz, Lx=LX, Lz=Lz, dt=dt, K=K, eps=pp['eps'], m=pp['m'],
                           c=pp['c'], regime=pp['regime'], four_eps2_m=pp['four_eps2_m'], dx=dx, dz=dz,
                           gamma=pp['c'] / dx ** 2, eta=pp['m'] * dz ** 2 / dx ** 2, patch=kind,
                           patch_size=(14 * Nz - 4 if kind == 'star' else 10 * Nz - 2),
                           its_v=iv, conv_v=int(cv), its_full=ifu, conv_full=int(cf))
                for k in ['ndof', 'lmin', 'lmax', 'cond', 'n_small', 'frac_small', 'lmin_nh', 'cond_nh',
                          'lmin_hf', 'cond_hf', 'n_small_hf', 'frac_small_hf']:
                    row[k] = s[k]
                row['seconds'] = time.time() - t0
                w.writerow({k: (f'{v:.10g}' if isinstance(v, float) else v) for k, v in row.items()})
                f.flush()
    print(f'{fname} written in {time.time() - t00:.1f} s')


if __name__ == '__main__':
    run(sys.argv[1])
