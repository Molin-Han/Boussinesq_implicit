"""
modes.py -- spectra of the additive Schwarz operators through the exact vertical
mode decomposition of Section 5 of theory_kernel.md.

For every vertical mode sigma in {0, sigma_1, ..., sigma_{k nz - 1}} the 2D
operator and all patch subspaces (which consist of whole columns) reduce to the
1D problem

    a_sigma((U,W),(U,W)) = |U|^2 + m |W|^2 + c |U' + sigma W|^2,
    U in CG_k(x) periodic,  W in DG_{k-1}(x)       (no W for sigma = 0).

spectrum(2D preconditioned operator) = union over the modes of the spectra of
the 1D preconditioned operators.
"""
import numpy as np
import scipy.linalg as sla
from fe1d import (Space1D, vertical_modes, mode_operator, mode_patches,
                  mode_prolongation, asm_inverse, gen_eigs)


def parameters(Lx, Lz, dt, N, K):
    """m, c, regime as in the experiments: delta = K * dt / (2 C_M^2), c = dt/(2 delta) = C_M^2 / K"""
    m = 1.0 + (dt * N) ** 2 / 4.0
    CM2 = max(Lx ** 2 / (4 * np.pi ** 2), m * Lz ** 2 / np.pi ** 2)
    c = CM2 / K
    eps = Lz / Lx
    regime = "weak" if 4 * eps ** 2 * m <= 1 else "strong"
    return m, c, regime


def level_spaces(Nx_fine, Lx, k, nlevels):
    """x-spaces, coarsest first; Nx halves from level to level"""
    out = []
    for l in range(nlevels):
        Nx = Nx_fine // 2 ** (nlevels - 1 - l)
        out.append(Space1D(Nx, Lx / Nx, k, periodic=True))
    return out


def additive_spectra(k, Nx, nz, Lx, Lz, m, c, kind, nlevels=3, sigmas=None, return_all=False):
    """For each vertical mode: eigenvalues of
         one-level   : ASM_fine A
         two-level   : (ASM_fine + P A_H^{-1} P^T) A            (H = 2h)
         multilevel  : (sum_{l>=1} P_l ASM_l P_l^T + P_0 A_0^{-1} P_0^T) A   (nlevels levels)
    Returns dict with arrays over the modes (mode 0 = sigma 0 first)."""
    dz = Lz / nz
    if sigmas is None:
        sigmas = np.concatenate([[0.0], vertical_modes(nz, dz, k)])
    X = level_spaces(Nx, Lx, k, nlevels)
    res = {"sigma": np.array(sigmas), "one": [], "two": [], "multi": []}
    for s in sigmas:
        A = mode_operator(X[-1], s, m, c)
        Minv = asm_inverse(A, mode_patches(X[-1], kind, s))
        e1 = gen_eigs(Minv, A)
        # two level
        P = mode_prolongation(X[-2], X[-1], s)
        AH = P.T @ A @ P
        B2 = Minv + P @ np.linalg.solve(AH, P.T)
        e2 = gen_eigs(B2, A)
        # multilevel (Galerkin coarse operators = rediscretised ones, spaces are nested)
        Bm = Minv.copy()
        Pl = np.eye(A.shape[0])
        for l in range(nlevels - 2, -1, -1):
            Pl = Pl @ mode_prolongation(X[l], X[l + 1], s)
            Al = Pl.T @ A @ Pl
            if l == 0:
                Bm += Pl @ np.linalg.solve(Al, Pl.T)
            else:
                Bm += Pl @ asm_inverse(Al, mode_patches(X[l], kind, s)) @ Pl.T
        em = gen_eigs(Bm, A)
        res["one"].append(e1)
        res["two"].append(e2)
        res["multi"].append(em)
    if return_all:
        return res
    out = {"sigma": res["sigma"]}
    for key in ("one", "two", "multi"):
        out[key + "_min"] = np.array([e[0] for e in res[key]])
        out[key + "_max"] = np.array([e[-1] for e in res[key]])
    return out


def summarise(out):
    """(lambda_min over baroclinic modes n>=1, lambda_min of barotropic mode n=0, lambda_max overall)
    for one/two/multi"""
    s = {}
    for key in ("one", "two", "multi"):
        lm = out[key + "_min"]
        s[key] = dict(min_bt=lm[0], min_bc=lm[1:].min(), argmin_bc=1 + int(np.argmin(lm[1:])),
                      min=lm.min(), max=out[key + "_max"].max())
    return s
