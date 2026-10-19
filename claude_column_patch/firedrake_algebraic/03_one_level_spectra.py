'''
03  One-level spectra of the additive Schwarz preconditioned velocity operator.

Operator (checked against the solver objects by 01_verify_operators_against_solver.py)

    S = (u_x, w_x) + m (u_z, w_z) + c (div u, div w),  m = 1 + (dt N)^2/4,  c = dt/(2 delta),

on the free velocity DoFs, delta = K delta_thr(dt, L_x, L_z) exactly as in the sweeps.
Preconditioner  P^{-1} = sum_j R_j^T (R_j S R_j^T)^{-1} R_j  with the velocity DoFs of the
patches of firedrake.ASMStarPC (construct_dim 0) and firedrake.ASMVankaPC (construct_dim 1).
(05_full_vs_reduced.py shows that the spectrum of the ASM-preconditioned FULL field-1
operator on (u, u_y, b) is this spectrum plus the eigenvalue 1 (vanka) or 2 (star).)

For every case: extreme eigenvalues, condition number, number of eigenvalues below
0.5, 0.1, 0.01, 0.001, and for the eigenvectors with eigenvalue < 0.1 the mean fraction of
their energy that is NOT divergence (1 means divergence free).

usage: python 03_one_level_spectra.py SWEEP [--degree 2]
   SWEEP in  K | ar | ncol | nz | dt | window
output: results_03_<SWEEP>_deg<k>.csv
'''
import os
import sys
import time
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fa_common as fa
from fa_common import *          # noqa

parser = argparse.ArgumentParser()
parser.add_argument("sweep", choices=["K", "ar", "ncol", "nz", "dt", "window"])
parser.add_argument("--degree", type=int, default=2)
args, _ = parser.parse_known_args()
deg = args.degree

L0, H0 = fa.NOMINAL_LENGTH, fa.NOMINAL_HEIGHT
KS = [0.1, 1.0, 10.0, 100.0, 1000.0]
DTS = [1., 10., 50., 100., 500., 1e3, 5e3, 1e4, 5e4, 1e5]
HEIGHTS = [12e3, 8e3, 4e3, 2e3, 1e3, 5e2, 2.5e2]

header = ["sweep", "degree", "ncol", "nz", "L_x", "L_z", "dx", "dz", "dt", "K", "delta", "m", "c",
          "c_over_dx2", "c_over_dz2", "regime", "patch", "n", "lambda_min", "lambda_max", "cond",
          "n_lt_0.5", "n_lt_0.1", "n_lt_0.01", "n_lt_0.001",
          "divfree_fraction_of_modes_lt_0.1", "seconds"]
rows = []
_cache = {}


def problem(nx, nz, length, height):
    key = (nx, nz, length, height, deg)
    if key not in _cache:
        _cache.clear()          # keep the memory footprint small
        _cache[key] = fa.VelocityProblem(nx, nz, length, height, deg)
    return _cache[key]


def run_case(sweep, nx, nz, length, height, dt, K, shift_geometry=None):
    '''shift_geometry = (L_x, L_z) used in the formula of delta (default: the mesh geometry).'''
    pb = problem(nx, nz, length, height)
    Ls, Hs = (length, height) if shift_geometry is None else shift_geometry
    delta, m, c = pb.params(dt, K, Ls, Hs)
    S = pb.S(m, c)
    Bd = pb.B.toarray()
    Mpinv = np.linalg.inv(pb.Mp.toarray())
    Dmat = c * (Bd.T @ Mpinv @ Bd)      # = c (div u, div w) since div V is contained in Q
    regime = cps.regime(dt, Ls, Hs)
    for kind in ["star", "vanka"]:
        t0 = time.time()
        Pinv = fa.asm_inverse_dense(S, pb.patches[kind])
        lam, X = fa.spd_preconditioned_eigs(S, Pinv, return_vectors=True)
        s = fa.spectrum_summary(lam)
        small = lam < 0.1
        if small.any():
            Xs = X[:, small]
            # X is S-normalised: x^T S x = 1 ; divergence part of the energy = x^T Dmat x
            divpart = np.einsum("ij,ij->j", Xs, Dmat @ Xs)
            divfree = float(np.mean(1.0 - divpart))
        else:
            divfree = float("nan")
        sec = time.time() - t0
        rows.append([sweep, deg, pb.ncol, nz, length, height, pb.dx, pb.dz, dt, K, delta, m, c,
                     c / pb.dx ** 2, c / pb.dz ** 2, regime, kind, s["n"], s["lmin"], s["lmax"], s["cond"],
                     s["n_lt_0.5"], s["n_lt_0.1"], s["n_lt_0.01"], s["n_lt_0.001"], divfree, sec])
        PETSc.Sys.Print(f"[{sweep}] deg {deg} {pb.ncol:3d}x{nz:<3d} L_z={height:7.0f} dt={dt:8.0f} K={K:7.1f} "
                        f"c/dx^2={c / pb.dx ** 2:9.3e} c/dz^2={c / pb.dz ** 2:9.3e} {regime:6s} {kind:5s} n={s['n']:5d} "
                        f"lmin={s['lmin']:.3e} lmax={s['lmax']:.3f} cond={s['cond']:.3e} "
                        f"#<0.5={s['n_lt_0.5']:5d} #<0.1={s['n_lt_0.1']:5d} #<0.01={s['n_lt_0.01']:5d} "
                        f"#<0.001={s['n_lt_0.001']:5d} divfree(<0.1)={divfree:.3f} ({sec:.0f}s)")


if args.sweep == "K":
    # nominal geometry, dt = 50 s, K sweep, on three small meshes
    for (nx, nz) in [(4, 8), (8, 8), (4, 16)]:
        for K in KS:
            run_case("K", nx, nz, L0, H0, 50.0, K)
elif args.sweep == "ar":
    # aspect ratio sweep at K = 1: L_z varies, delta = delta_thr of each geometry
    for dt in [50.0, 1.0e4]:
        for H in HEIGHTS:
            run_case("ar", 4, 8, L0, H, dt, 1.0)
elif args.sweep == "ncol":
    for nx in [2, 4, 8, 16]:
        for K in [1.0, 100.0]:
            run_case("ncol", nx, 6, L0, H0, 50.0, K)
elif args.sweep == "nz":
    for nz in [4, 8, 16, 32]:
        for dt in [50.0, 1.0e5]:
            run_case("nz", 4, nz, L0, H0, dt, 1.0)
elif args.sweep == "dt":
    for dt in DTS:
        for K in [1.0, 100.0]:
            run_case("dt", 4, 8, L0, H0, dt, K)
elif args.sweep == "window":
    # a window of the NOMINAL mesh: same cells (dx = 3750 m, dz = 200 m, 20 layers), same
    # delta as the nominal 300 km x 4 km problem, but only 16 or 32 columns (periodic)
    dxn = L0 / 80.0
    for nx in ([4, 8] if deg == 2 else [4, 8, 16]):
        ncol = 4 * nx
        for K in KS:
            run_case("window", nx, 20, ncol * dxn, H0, 50.0, K, shift_geometry=(L0, H0))

fa.write_table(os.path.join(HERE, f"results_03_{args.sweep}_deg{deg}.csv"), header, rows)
