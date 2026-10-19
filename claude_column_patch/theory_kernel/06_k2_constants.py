"""
06_k2_constants.py

Degree k = 2 (the degree of the solver).  Constants of Section 7 of theory_kernel.md.

PART 1  C_inf(rho) = sup_U [ S_rho(I_H U) + S_rho(U - I_H U) ] / S_rho(U),
        S_rho(V) = |V|^2 + rho h^2 |V'|^2,  U in CG_2(h) periodic,
        I_H U in CG_2(2h) the interpolant of U at the vertices of the FINE mesh
        (these are exactly the nodes of CG_2(2h)), so that U - I_H U is a sum of
        fine cell bubbles.  This is the two-level stable-decomposition constant of
        the single-column patch in the limit c -> infinity (Theorem 7.2).

PART 2  two-level lambda_min of the single-column patch for growing gamma = c/dx^2
        and several rho: converges to 1 / C_inf(rho).

PART 3  one-level single-column patch: K_0 / gamma for growing gamma (thin limit: 15).

PART 4  number of levels: lambda_min of the multilevel additive operator for
        2, 3, 4, 5 levels (coarsest level solved exactly), star and column, k = 1, 2.

Run: python 06_k2_constants.py > 06_k2_constants.txt
"""
import numpy as np
import scipy.linalg as sla
from fe1d import Space1D, prolongation_x, mode_operator, mode_patches, mode_prolongation, asm_inverse, gen_eigs


def interp_fine_vertices(XH, Xh):
    """I_H : CG_2(h) -> CG_2(2h) -> CG_2(h), interpolation at the fine vertices"""
    N = Xh.n
    Pc, _ = prolongation_x(XH, Xh)
    R = np.zeros((XH.nc, Xh.nc))
    for J in range(N // 2):
        R[J, 2 * J] = 1.0                    # coarse vertex J   = fine vertex 2J
        R[XH.nvert + J, 2 * J + 1] = 1.0     # coarse midpoint J = fine vertex 2J+1
    return Pc @ R


def C_inf(rho, N=16):
    h = 1.0
    Xh, XH = Space1D(N, h, 2, True), Space1D(N // 2, 2 * h, 2, True)
    IH = interp_fine_vertices(XH, Xh)
    Id = np.eye(Xh.nc)
    if np.isinf(rho):
        # H1 seminorm only: work on the complement of the constants
        S = Xh.Kc + 1e-12 * Xh.Mc
    else:
        S = Xh.Mc + rho * h * h * Xh.Kc
    T = IH.T @ S @ IH + (Id - IH).T @ S @ (Id - IH)
    ev = sla.eigh(T, S, eigvals_only=True)
    return ev[-1]


def two_level_min(kind, k, N, g, tau, m, nlev=2):
    h = 1.0
    X = [Space1D(N // 2 ** (nlev - 1 - l), h * 2 ** (nlev - 1 - l), k, True) for l in range(nlev)]
    out = []
    for sig in ((0.0,) if tau is None else (tau / h,)):
        A = mode_operator(X[-1], sig, m, g * h * h)
        B = asm_inverse(A, mode_patches(X[-1], kind, sig))
        Pl = np.eye(A.shape[0])
        for l in range(nlev - 2, -1, -1):
            Pl = Pl @ mode_prolongation(X[l], X[l + 1], sig)
            Al = Pl.T @ A @ Pl
            if l == 0:
                B += Pl @ np.linalg.solve(Al, Pl.T)
            else:
                B += Pl @ asm_inverse(Al, mode_patches(X[l], kind, sig)) @ Pl.T
        e = gen_eigs(B, A)
        out.append((e[0], e[-1]))
    return out[0]


def main():
    print("PART 1: C_inf(rho) for k = 2 (independent of Nx: local constant)")
    for N in (8, 16, 32):
        print(f"  Nx = {N:3d}: " + ", ".join(f"rho={r:g}: {C_inf(r, N):.4f}" for r in (0.0, 1e-3, 1e-2, 0.1, 1.0, 10.0, 1e3, np.inf)))

    print("\nPART 2: single-column patch, k = 2, two level: lambda_min for growing gamma (tau fixed, m = 1)")
    print("        rho = gamma m/(m + gamma tau^2) -> m/tau^2 ;  limit 1/C_inf(m/tau^2)")
    for tau in (60.0, 3.0, 1.0, 0.3):
        m = 1.0
        row = []
        for g in (1.0, 10.0, 1e2, 1e3, 1e4, 1e6):
            lmin, lmax = two_level_min("vanka", 2, 16, g, tau, m)
            row.append(f"{lmin:.4f}")
        print(f"  tau = {tau:5.1f} (rho_inf = {m/tau**2:8.3g}): gamma = 1, 10, 1e2, 1e3, 1e4, 1e6 -> "
              + " ".join(row) + f"   | 1/C_inf = {1/C_inf(m/tau**2):.4f}")
    row = []
    for g in (1.0, 10.0, 1e2, 1e3, 1e4, 1e6):
        lmin, lmax = two_level_min("vanka", 2, 16, g, None, 1.0)
        row.append(f"{lmin:.4f}")
    print("  barotropic mode (sigma = 0, rho = gamma)   : gamma = 1, 10, 1e2, 1e3, 1e4, 1e6 -> " + " ".join(row))
    print("  same for the star patch:")
    for tau in (60.0, 3.0, 1.0, 0.3):
        row = []
        for g in (1.0, 10.0, 1e2, 1e3, 1e4, 1e6):
            lmin, lmax = two_level_min("star", 2, 16, g, tau, 1.0)
            row.append(f"{lmin:.4f}")
        print(f"  tau = {tau:5.1f}: " + " ".join(row))
    row = []
    for g in (1.0, 10.0, 1e2, 1e3, 1e4, 1e6):
        lmin, lmax = two_level_min("star", 2, 16, g, None, 1.0)
        row.append(f"{lmin:.4f}")
    print("  barotropic    : " + " ".join(row))

    print("\nPART 3: single-column patch, one level: K_0 / gamma = 1 / (gamma lambda_min), tau = 60, m = 1 (thin limit)")
    for k in (1, 2):
        row = []
        for g in (10.0, 1e2, 1e3, 1e4):
            h = 1.0
            X = Space1D(32, h, k, True)
            A = mode_operator(X, 60.0, 1.0, g)
            e = gen_eigs(asm_inverse(A, mode_patches(X, "vanka", 60.0)), A)
            row.append(f"{1/(g*e[0]):.3f}")
        print(f"  k = {k}: gamma = 10, 1e2, 1e3, 1e4 -> K_0/gamma = " + " ".join(row)
              + ("   (theory: 9/2)" if k == 1 else "   (theory: 15)"))

    print("\nPART 4: number of levels (finest Nx = 64, gamma = 1e4 on the finest level, tau = 60, m = 1):"
          " lambda_min / lambda_max of the multilevel additive operator")
    for k in (1, 2):
        for kind in ("star", "vanka"):
            row = []
            for nlev in (2, 3, 4, 5):
                lmin, lmax = two_level_min(kind, k, 64, 1e4, 60.0, 1.0, nlev=nlev)
                lmin0, lmax0 = two_level_min(kind, k, 64, 1e4, None, 1.0, nlev=nlev)
                row.append(f"{nlev} levels: bc {lmin:.3e}/{lmax:.2f} bt {lmin0:.3e}/{lmax0:.2f}")
            print(f"  k = {k} {kind:5s}: " + " | ".join(row))


if __name__ == "__main__":
    main()
