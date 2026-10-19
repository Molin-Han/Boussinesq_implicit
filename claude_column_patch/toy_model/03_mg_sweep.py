'''
Iteration counts of FGMRES (relative tolerance 1e-8, random right-hand side, restart 30)
preconditioned by ONE multigrid cycle for the toy H(div) problem.

Multigrid: 3 levels (as --refinement 2 in the solver script), coarsening in x only by
factors of 2, the same N_z on every level, sparse LU on the coarsest level, smoother =
6 iterations of GMRES preconditioned by the additive Schwarz method (star or single-column
patch), pre and post smoothing.  N_x is the number of columns on the FINEST level.

usage:  python 03_mg_sweep.py <group> <cycle>
    group  geom : N_x = 80, N_z = 20,  L_z x dt x K       (aspect ratio, time step, shift)
           nx   : N_z = 20, L_z = 4 km, N_x x dt x K
           nz   : N_x = 80, L_z = 4 km, N_z x dt x K
    cycle  v    : V cycle (the task specification)
           full : PETSc -pc_mg_type full with V cycles (what the solver script uses)

Output: tables/mg_<group>_<cycle>.csv
'''
import csv
import os
import sys
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
NOMINAL = dict(Nx=80, Nz=20, Lz=4e3)

FIELDS = ['group', 'cycle', 'Nx', 'Nz', 'Lx', 'Lz', 'dt', 'K', 'eps', 'm', 'c', 'regime', 'four_eps2_m',
          'dx', 'dz', 'gamma', 'eta', 'patch', 'nlevels', 'smooth_its', 'rtol', 'seed',
          'its', 'converged', 'true_rel_res', 'mean_reduction_per_it', 'seconds']


def configs(group):
    if group == 'geom':
        return [(NOMINAL['Nx'], NOMINAL['Nz'], Lz, dt, K) for Lz, dt, K in itertools.product(LZS, DTS, KS)]
    if group == 'nx':
        return [(Nx, NOMINAL['Nz'], NOMINAL['Lz'], dt, K) for Nx, dt, K in itertools.product(NXS, DTS, KS)]
    if group == 'nz':
        return [(NOMINAL['Nx'], Nz, NOMINAL['Lz'], dt, K) for Nz, dt, K in itertools.product(NZS, DTS, KS)]
    raise ValueError(group)


def run(group, cycle, nlevels=3, smooth_its=6, rtol=1e-8, seed=0, maxit=300):
    os.makedirs('tables', exist_ok=True)
    fname = f'tables/mg_{group}_{cycle}.csv'
    t00 = time.time()
    with open(fname, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for (Nx, Nz, Lz, dt, K) in configs(group):
            pp = rt.physical_parameters(LX, Lz, dt, K, N=NFREQ)
            dx, dz = LX / Nx, Lz / Nz
            for kind in rt.PATCHES:
                t0 = time.time()
                its, conv, res, hist = rt.mg_fgmres_iterations(
                    Nx, Nz, LX, Lz, pp['m'], pp['c'], kind, nlevels=nlevels, smooth_its=smooth_its,
                    cycle=cycle, rtol=rtol, seed=seed, maxit=maxit, return_history=True)
                red = (hist[-1] / hist[0]) ** (1.0 / max(its, 1))
                row = dict(group=group, cycle=cycle, Nx=Nx, Nz=Nz, Lx=LX, Lz=Lz, dt=dt, K=K, eps=pp['eps'],
                           m=pp['m'], c=pp['c'], regime=pp['regime'], four_eps2_m=pp['four_eps2_m'],
                           dx=dx, dz=dz, gamma=pp['c'] / dx ** 2, eta=pp['m'] * dz ** 2 / dx ** 2,
                           patch=kind, nlevels=nlevels, smooth_its=smooth_its, rtol=rtol, seed=seed,
                           its=its, converged=int(conv), true_rel_res=res, mean_reduction_per_it=red,
                           seconds=time.time() - t0)
                w.writerow({k: (f'{v:.10g}' if isinstance(v, float) else v) for k, v in row.items()})
                f.flush()
    print(f'{fname} written in {time.time() - t00:.1f} s')


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])
