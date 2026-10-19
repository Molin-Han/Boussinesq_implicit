"""
04_spectra_parameter_study.py

Smallest / largest eigenvalues of the additive Schwarz preconditioned velocity
operator for the star patch and the single-column (vanka) patch, degree k = 1, 2,
computed exactly through the vertical mode decomposition (modes.py), for the
sweeps of the experiments:

    K  in {0.01, 0.1, 1, 10, 100, 1000}           (c = C_M^2 / K)
    Lz in {12, 8, 4, 2, 1, 0.5, 0.25} km
    nx in {10, 20, 40, 60, 80}  (Nx = 4 nx columns on the finest level)
    nz in {10, 20, 40, 60, 80}
    dt in {1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5} s

Reported for every case
    one-level   lambda_min over the baroclinic modes (n >= 1) and of the barotropic mode (n = 0)
    two-level   lambda_min (coarse space: coarsening by 2 in x, exact coarse solve)
    three-level lambda_min (the hierarchy of the solver: refinement = 2, exact solve on the coarsest level)
    lambda_max  of the three-level operator
The stable-decomposition constant is K_0 = 1 / lambda_min.

Writes 04_spectra_parameter_study.txt and 04_spectra_parameter_study.csv
Run: python 04_spectra_parameter_study.py
"""
import sys
import numpy as np
from modes import parameters, additive_spectra, summarise

Lx0, Lz0, N0, nx0, nz0, dt0, K0 = 300e3, 4e3, 1e-2, 20, 20, 50.0, 1.0
REF = 2

rows = []
out = open("04_spectra_parameter_study.txt", "w")


def P(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    out.write(s + "\n")
    out.flush()


def case(label, value, k, Lx=Lx0, Lz=Lz0, nx=nx0, nz=nz0, dt=dt0, K=K0):
    Nx = nx * 2 ** REF
    m, c, regime = parameters(Lx, Lz, dt, N0, K)
    dx, dz = Lx / Nx, Lz / nz
    line = f"{label:>4s}={value:<8g} k={k} g=c/dx^2={c/dx**2:9.3e} c/dz^2={c/dz**2:9.3e} m={m:9.3e} {regime:6s}|"
    for kind in ("star", "vanka"):
        s = summarise(additive_spectra(k, Nx, nz, Lx, Lz, m, c, kind, nlevels=REF + 1))
        line += (f" {kind}: 1L bc {s['one']['min_bc']:.2e} bt {s['one']['min_bt']:.2e}"
                 f" 2L {s['two']['min']:.2e} 3L {s['multi']['min']:.2e} max {s['multi']['max']:.2f} |")
        rows.append(dict(sweep=label, value=value, k=k, patch=kind, Nx=Nx, nz=nz, Lx=Lx, Lz=Lz, dt=dt, K=K,
                         m=m, c=c, gamma=c / dx ** 2, c_over_dz2=c / dz ** 2, regime=regime,
                         one_min_bc=s['one']['min_bc'], one_min_bt=s['one']['min_bt'],
                         one_max=s['one']['max'],
                         two_min=s['two']['min'], two_max=s['two']['max'],
                         three_min=s['multi']['min'], three_max=s['multi']['max'],
                         two_argmin_mode=int(np.argmin(np.r_[s['two']['min_bt'], s['two']['min_bc']]))))
    P(line)


P("lambda_min of the additive Schwarz preconditioned operator; K_0 = 1/lambda_min")
P("1L = one level (bc: min over baroclinic modes n>=1, bt: barotropic mode n=0), 2L = two level, 3L = three level")
P("nominal: Lx=300 km, Lz=4 km, N=1e-2, nx=20 (Nx=80), nz=20, dt=50, K=1\n")

which = sys.argv[1:] or ["K", "Lz", "nx", "nz", "dt"]
for k in (1, 2):
    if "K" in which:
        P(f"--- sweep K (shift factor), degree {k}")
        for K in (0.01, 0.1, 1, 10, 100, 1000):
            case("K", K, k, K=K)
    if "Lz" in which:
        P(f"--- sweep Lz [km] (aspect ratio), degree {k}")
        for Lz in (12, 8, 4, 2, 1, 0.5, 0.25):
            case("Lz", Lz, k, Lz=Lz * 1e3)
    if "nx" in which:
        P(f"--- sweep nx (columns on the finest level = 4 nx), degree {k}")
        for nx in (10, 20, 40, 60, 80):
            case("nx", nx, k, nx=nx)
    if "nz" in which:
        P(f"--- sweep nz, degree {k}")
        for nz in (10, 20, 40, 60, 80):
            case("nz", nz, k, nz=nz)
    if "dt" in which:
        P(f"--- sweep dt [s], degree {k}")
        for dt in (1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5):
            case("dt", dt, k, dt=dt)

import csv
with open("04_spectra_parameter_study.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
out.close()
