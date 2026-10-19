"""
03_mode_decoupling_check.py

Checks Lemma 5.1 of theory_kernel.md (exact decoupling into vertical modes):
the spectrum of the additive Schwarz preconditioned 2D operator equals the union
over the vertical modes sigma_n of the spectra of the 1D mode problems.

Run: python 03_mode_decoupling_check.py > 03_mode_decoupling_check.txt
"""
import numpy as np
from fe1d import Model2D, asm_inverse, gen_eigs
from modes import additive_spectra

Lx, Lz = 300e3, 4e3
ok = True
for k in (1, 2):
    for (Nx, nz, m, c) in ((8, 3, 1.0625, 5.0e8), (8, 4, 40.0, 2.3e9), (16, 3, 1.0625, 2.0e7)):
        fine = Model2D(Nx, nz, Lx, Lz, k)
        coarse = Model2D(Nx // 2, nz, Lx, Lz, k)
        A = fine.operator(m, c)
        P = fine.prolongation_from(coarse)
        AH = P.T @ A @ P
        for kind in ("star", "vanka"):
            Minv = asm_inverse(A, fine.patches(kind))
            e1 = gen_eigs(Minv, A)
            e2 = gen_eigs(Minv + P @ np.linalg.solve(AH, P.T), A)
            r = additive_spectra(k, Nx, nz, Lx, Lz, m, c, kind, nlevels=2, return_all=True)
            m1 = np.sort(np.concatenate(r["one"]))
            m2 = np.sort(np.concatenate(r["two"]))
            d1 = np.max(np.abs(m1 - e1) / np.abs(e1))
            d2 = np.max(np.abs(m2 - e2) / np.abs(e2))
            ok &= (len(m1) == len(e1)) and d1 < 1e-7 and d2 < 1e-7
            print(f"k={k} Nx={Nx} nz={nz} m={m:g} c/dx^2={c/(Lx/Nx)**2:.3g} {kind:5s}: "
                  f"#eig 2D = {len(e1)}, #eig modes = {len(m1)}, "
                  f"max rel diff one-level {d1:.1e}, two-level {d2:.1e}")
print("union of the mode spectra = spectrum of the 2D operator:", ok)
