"""
05_lfa_k1_closed_forms.py

Degree k = 1.  Checks the closed formulas of Section 6 of theory_kernel.md for the
spectrum of the one-level additive Schwarz operator of a vertical mode sigma
against the eigenvalues computed from the matrices, and the rigorous two-level
bounds for the single-column patch.

Notation:  g = gamma = c/dx^2,  tau = sigma dx,  rho = g m / (m + g tau^2) = c_eff/dx^2,
           theta = 2 pi j / Nx,  C = cos(theta/2), s = sin(theta/2), mu = (2 + cos theta)/3.

star, baroclinic mode   : eigenvalues {2} and (mu + 4 rho s^2) / (2/3 + 2 rho)
star, barotropic mode   : (mu + 4 g s^2) / (2/3 + 2 g)
column, baroclinic mode : the two eigenvalues of a Hermitian 2 x 2 matrix with
       det = shat [ 2 (a+ C^2 + a- s^2) + 4 G a+ a- C^2 ]
       tr  = 1 + 2 shat (a+ C^2 + a- s^2) + 2 G (a+ sin^2 theta + a- cos^2 theta)
       shat = mu + 4 rho s^2, G = g - rho, a+ = 1/(5/6 + g), a- = 1/(1/2 + g + 2 rho)
column, barotropic mode : (mu + 4 g s^2) 2 (C^2/(5/6 + g) + s^2/(1/2 + 3 g))

Run: python 05_lfa_k1_closed_forms.py > 05_lfa_k1_closed_forms.txt
"""
import numpy as np
from fe1d import Space1D, mode_operator, mode_patches, mode_prolongation, asm_inverse, gen_eigs


def closed_forms(N, g, tau, m):
    th = 2 * np.pi * np.arange(N) / N
    C, s = np.cos(th / 2), np.sin(th / 2)
    mu = (2 + np.cos(th)) / 3
    out = {}
    out[("star", "bt")] = np.sort((mu + 4 * g * s ** 2) / (2 / 3 + 2 * g))
    out[("vanka", "bt")] = np.sort((mu + 4 * g * s ** 2) * 2 * (C ** 2 / (5 / 6 + g) + s ** 2 / (0.5 + 3 * g)))
    rho = g * m / (m + g * tau ** 2)
    shat = mu + 4 * rho * s ** 2
    out[("star", "bc")] = np.sort(np.concatenate([2 * np.ones(N), shat / (2 / 3 + 2 * rho)]))
    G = g - rho
    ap, am = 1 / (5 / 6 + g), 1 / (0.5 + g + 2 * rho)
    det = shat * (2 * (ap * C ** 2 + am * s ** 2) + 4 * G * ap * am * C ** 2)
    tr = 1 + 2 * shat * (ap * C ** 2 + am * s ** 2) + 2 * G * (ap * np.sin(th) ** 2 + am * np.cos(th) ** 2)
    disc = np.sqrt(np.maximum(tr ** 2 - 4 * det, 0))
    lp = 0.5 * (tr + disc)
    lm = det / lp
    out[("vanka", "bc")] = np.sort(np.concatenate([lm, lp]))
    return out, rho


def main():
    print("PART 1: closed formulas (k = 1, one level) versus matrices")
    print(f"{'Nx':>4s} {'gamma':>10s} {'tau':>9s} {'m':>9s} {'rho':>10s} | max relative deviation: "
          f"star-bt  star-bc  column-bt  column-bc")
    ok = True
    for (N, g, tau, m) in ((16, 1.6, 5.9, 1.0625), (32, 162.0, 2.9, 1.0625), (32, 162.0, 65.0, 1.0625),
                           (40, 2594.0, 0.74, 1.0625), (32, 28820.0, 65.0, 250001.0), (32, 0.16, 12.0, 7.25),
                           (24, 40.0, 1e-3, 1.0)):
        h = 1.0
        X = Space1D(N, h, 1, True)
        cf, rho = closed_forms(N, g, tau, m)
        devs = []
        for kind in ("star", "vanka"):
            for tag, sig in (("bt", 0.0), ("bc", tau / h)):
                A = mode_operator(X, sig, m, g * h * h)
                e = gen_eigs(asm_inverse(A, mode_patches(X, kind, sig)), A)
                d = np.max(np.abs(e - cf[(kind, tag)]) / np.abs(cf[(kind, tag)]))
                devs.append(d)
        ok &= max(devs) < 1e-6
        print(f"{N:4d} {g:10.4g} {tau:9.4g} {m:9.4g} {rho:10.4g} |   "
              + "  ".join(f"{d:8.1e}" for d in devs))
    print("closed formulas reproduce the computed spectra:", ok)

    print("\nPART 2: two-level operator, single-column patch, k = 1 (coarse space: coarsening by 2 in x)")
    print("  lam_-(pi)  = smaller root of  x^2 - tr x + det  at theta = pi:  det = 2 shat a-,  tr = 1 + 2 a- (shat + G),"
          " shat = 1/3 + 4 rho")
    print("  rigorous bounds of Theorem 6.3:   (3 g / 2) / (1 + 12 m / tau^2)  <=  K_0^TL  <=  2 + 7.5 g")
    print(f"{'Nx':>4s} {'gamma':>10s} {'tau':>9s} {'m':>9s} | {'lam_min 1L':>11s} {'lam_min 2L':>11s} {'lam_-(pi)':>11s}"
          f" {'2(1+12rho)/(9g)':>16s} | {'lower bd':>9s} {'K_0^TL':>9s} {'upper bd':>9s}")
    ok2 = True
    for (N, g, tau, m) in ((16, 1.6, 5.9, 1.0625), (32, 16.2, 2.9, 1.0625), (32, 162.0, 2.9, 1.0625),
                           (32, 162.0, 65.0, 1.0625), (64, 2594.0, 65.0, 1.0625), (32, 28820.0, 65.0, 250001.0),
                           (32, 1621.0, 65.0, 626.0)):
        h = 1.0
        Xh, XH = Space1D(N, h, 1, True), Space1D(N // 2, 2 * h, 1, True)
        sig = tau / h
        A = mode_operator(Xh, sig, m, g * h * h)
        Minv = asm_inverse(A, mode_patches(Xh, "vanka", sig))
        P = mode_prolongation(XH, Xh, sig)
        e1 = gen_eigs(Minv, A)
        e2 = gen_eigs(Minv + P @ np.linalg.solve(P.T @ A @ P, P.T), A)
        rho = g * m / (m + g * tau ** 2)
        G = g - rho
        am = 1 / (0.5 + g + 2 * rho)
        shat = 1 / 3 + 4 * rho
        det, tr = 2 * shat * am, 1 + 2 * am * (shat + G)
        lpi = det / (0.5 * (tr + np.sqrt(tr ** 2 - 4 * det)))
        asym = 2 * (1 + 12 * rho) / (9 * g)
        lower = 1.5 * g / (1 + 12 * m / tau ** 2)
        upper = 2 + 7.5 * g
        K0 = 1 / e2[0]
        ok2 &= (lower <= K0 * (1 + 1e-9)) and (K0 <= upper * (1 + 1e-9)) and abs(e2[0] - lpi) / lpi < 1e-6
        print(f"{N:4d} {g:10.4g} {tau:9.4g} {m:9.4g} | {e1[0]:11.4e} {e2[0]:11.4e} {lpi:11.4e} {asym:16.4e} |"
              f" {lower:9.3g} {K0:9.3g} {upper:9.3g}")
    print("two-level lambda_min = lam_-(pi) and the bounds hold:", ok2)

    print("\nPART 3: two-level operator, star patch, k = 1: lambda_min (Theorem 6.2 proves >= 1/10)")
    for (N, g, tau, m) in ((32, 0.16, 2.9, 1.0625), (32, 162.0, 2.9, 1.0625), (32, 162.0, 65.0, 1.0625),
                           (64, 2594.0, 0.5, 1.0625), (32, 28820.0, 3.0, 250001.0), (64, 1e6, 0.05, 1.0)):
        h = 1.0
        Xh, XH = Space1D(N, h, 1, True), Space1D(N // 2, 2 * h, 1, True)
        vals = []
        for sig in (0.0, tau / h):
            A = mode_operator(Xh, sig, m, g * h * h)
            Minv = asm_inverse(A, mode_patches(Xh, "star", sig))
            P = mode_prolongation(XH, Xh, sig)
            e2 = gen_eigs(Minv + P @ np.linalg.solve(P.T @ A @ P, P.T), A)
            vals.append((e2[0], e2[-1]))
        rho = g * m / (m + g * tau ** 2)
        print(f"  Nx={N:3d} gamma={g:9.4g} tau={tau:6.3g} m={m:9.4g} rho={rho:9.3g}: barotropic lam in [{vals[0][0]:.4f}, {vals[0][1]:.4f}],"
              f" baroclinic lam in [{vals[1][0]:.4f}, {vals[1][1]:.4f}]")


if __name__ == "__main__":
    main()
