'''
Spectrum of the one-level additive Schwarz preconditioned operator B A for the star and the
single-column patch, full factorial sweep

    N_x in {20, 40, 80, 160, 320}, N_z in {10, 20, 40}, L_z in {12, 8, 4, 2, 1, 0.5, 0.25} km,
    dt in {1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5} s, K in {0.1, 1, 10, 100, 1000},
    L_x = 300 km, N = 1e-2 1/s, delta = K delta_thr.

All eigenvalues are computed from the 2 x 2 symbols (Fourier in x, sine/cosine in z), which
01_validate.py checks against the dense computation.

Output: tables/spectra_full_factorial.csv   (one row per configuration and patch)
'''
import csv
import os
import time
import itertools
import numpy as np
import rt0_model as rt

LX = 3.0e5
NFREQ = 1.0e-2
NXS = [20, 40, 80, 160, 320]
NZS = [10, 20, 40]
LZS = [12e3, 8e3, 4e3, 2e3, 1e3, 0.5e3, 0.25e3]
DTS = [1., 10., 50., 100., 500., 1e3, 5e3, 1e4, 5e4, 1e5]
KS = [0.1, 1., 10., 100., 1000.]

FIELDS = ['Nx', 'Nz', 'Lx', 'Lz', 'dt', 'K', 'N', 'eps', 'm', 'c', 'delta', 'regime', 'four_eps2_m',
          'dx', 'dz', 'gamma', 'eta', 'patch', 'patch_size', 'ndof',
          'lmin', 'lmax', 'cond', 'n_small', 'frac_small',
          'lam_uniform_flow', 'lmin_nh', 'cond_nh',
          'ndof_hf', 'lmin_hf', 'lmax_hf', 'cond_hf', 'n_small_hf', 'frac_small_hf',
          'argmin_theta_over_pi', 'argmin_l']


def row_for(Nx, Nz, Lz, dt, K, kind):
    pp = rt.physical_parameters(LX, Lz, dt, K, N=NFREQ)
    dx, dz = LX / Nx, Lz / Nz
    s = rt.spectral_summary(Nx, Nz, dx, dz, pp['m'], pp['c'], kind)
    row = dict(Nx=Nx, Nz=Nz, Lx=LX, Lz=Lz, dt=dt, K=K, N=NFREQ, eps=pp['eps'], m=pp['m'], c=pp['c'],
               delta=pp['delta'], regime=pp['regime'], four_eps2_m=pp['four_eps2_m'],
               dx=dx, dz=dz, gamma=pp['c'] / dx ** 2, eta=pp['m'] * dz ** 2 / dx ** 2,
               patch=kind, patch_size=(3 * Nz - 2 if kind == 'star' else 3 * Nz - 1))
    for k in ['ndof', 'lmin', 'lmax', 'cond', 'n_small', 'frac_small', 'lam_uniform_flow', 'lmin_nh',
              'cond_nh', 'ndof_hf', 'lmin_hf', 'lmax_hf', 'cond_hf', 'n_small_hf', 'frac_small_hf', 'argmin_l']:
        row[k] = s[k]
    row['argmin_theta_over_pi'] = s['argmin_theta'] / np.pi
    return row


if __name__ == '__main__':
    os.makedirs('tables', exist_ok=True)
    t0 = time.time()
    n = 0
    with open('tables/spectra_full_factorial.csv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for Nx, Nz, Lz, dt, K in itertools.product(NXS, NZS, LZS, DTS, KS):
            for kind in rt.PATCHES:
                w.writerow({k: (f'{v:.10g}' if isinstance(v, float) else v)
                            for k, v in row_for(Nx, Nz, Lz, dt, K, kind).items()})
                n += 1
    print(f'{n} rows written to tables/spectra_full_factorial.csv in {time.time() - t0:.1f} s')
