"""
07_mg_cycle_emulation.py

Emulation (mg_emulator.py) of the multigrid of the solver applied to the
symmetric velocity operator: number of FGMRES iterations (one full-multigrid
V(6,6) cycle with GMRES/additive-Schwarz smoothing per iteration) needed to
reduce the residual by 1e-6, for the star and the single-column patch and for
degree k = 1 and 2.

Two right-hand sides:  'phys' = load of the initial buoyancy perturbation of the
solver script,  'rand' = random vector (fixed seed).

usage: python 07_mg_cycle_emulation.py [K] [Lz] [nx] [nz] [dt] [ref]
       (no argument: all sweeps).  Output: 07_mg_cycle_emulation_<sweep>.txt
"""
import sys
import time
import numpy as np
from modes import parameters
from mg_emulator import Multigrid, fgmres, physical_rhs

Lx0, Lz0, N0, nx0, nz0, dt0, K0 = 300e3, 4e3, 1e-2, 20, 20, 50.0, 1.0


def run(k, kind, Lx=Lx0, Lz=Lz0, nx=nx0, nz=nz0, dt=dt0, K=K0, ref=2, maxit=150):
    m, c, regime = parameters(Lx, Lz, dt, N0, K)
    mg = Multigrid(k, nx, nz, Lx, Lz, m, c, kind, refinement=ref)
    mod = mg.levels[-1]["mod"]
    out = {}
    rng = np.random.default_rng(1)
    for name, b in (("phys", physical_rhs(mod)), ("rand", rng.standard_normal(mg.n))):
        x, hist = fgmres(mg.A, mg.fcycle, b, rtol=1e-6, maxit=maxit)
        its = len(hist) - 1
        conv = hist[-1] < 1e-6
        out[name] = (its if conv else f">{maxit}", hist[1], (hist[-1]) ** (1.0 / its))
    Nx = nx * 2 ** ref
    return out, c / (Lx / Nx) ** 2, m, regime


def sweep(label, values, fh, **fixed):
    def P(s):
        print(s)
        fh.write(s + "\n")
        fh.flush()
    P(f"\n=== sweep {label}  (other parameters nominal: Lx=300km Lz=4km nx=20 nz=20 dt=50 K=1 refinement=2)")
    P(f"{label:>6s} {'gamma':>10s} {'m':>10s} | " + " | ".join(
        f"k={k} {kind:5s} phys/rand" for k in (1, 2) for kind in ("star", "vanka")))
    for v in values:
        kw = dict(fixed)
        kw[label if label != "Lz" else "Lz"] = v * 1e3 if label == "Lz" else v
        cells = []
        for k in (1, 2):
            for kind in ("star", "vanka"):
                t0 = time.time()
                out, g, m, regime = run(k, kind, **kw)
                cells.append(f"{str(out['phys'][0]):>5s}/{str(out['rand'][0]):<5s} "
                             f"(r1 {out['phys'][1]:.1e})")
        P(f"{v:6g} {g:10.3e} {m:10.3e} | " + " | ".join(cells))


if __name__ == "__main__":
    which = sys.argv[1:] or ["K", "Lz", "nx", "nz", "dt", "ref"]
    for w in which:
        with open(f"07_mg_cycle_emulation_{w}.txt", "w") as fh:
            fh.write("FGMRES iterations to rtol 1e-6, one full-MG V(6,6) GMRES/ASM cycle per iteration; "
                     "'phys'/'rand' right-hand side; r1 = relative residual after the first iteration (phys)\n")
            if w == "K":
                sweep("K", (0.01, 0.1, 1, 10, 100, 1000), fh)
            elif w == "Lz":
                sweep("Lz", (12, 8, 4, 2, 1, 0.5, 0.25), fh)
            elif w == "nx":
                sweep("nx", (10, 20, 40, 60, 80), fh)
            elif w == "nz":
                sweep("nz", (10, 20, 40, 60, 80), fh)
            elif w == "dt":
                sweep("dt", (1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5), fh)
            elif w == "ref":
                # same finest mesh (Nx = 80), more levels
                def P(s):
                    print(s)
                    fh.write(s + "\n")
                    fh.flush()
                P("\n=== number of levels at fixed finest mesh Nx = 80 (nx_coarse = 80 / 2^ref)")
                P(f"{'ref':>6s} | " + " | ".join(f"k={k} {kind:5s} phys/rand" for k in (1, 2) for kind in ("star", "vanka")))
                for ref in (1, 2, 3, 4):
                    cells = []
                    for k in (1, 2):
                        for kind in ("star", "vanka"):
                            out, g, m, regime = run(k, kind, nx=80 // 2 ** ref, ref=ref)
                            cells.append(f"{str(out['phys'][0]):>5s}/{str(out['rand'][0]):<5s}")
                    P(f"{ref:6d} | " + " | ".join(cells))
