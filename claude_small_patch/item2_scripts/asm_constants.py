"""Throwaway: additive-Schwarz constants for the two 1D splittings behind decomposition (C).

(1) z-direction: CG_2^0(0,1), N_z cells, operator K + kappa^2 M (Dirichlet both ends),
    patches = 1D vertex stars {vertex l, bubble l-1/2, bubble l+1/2}, NO coarse space.
    This is what the 2D vertex stars impose on the hat-column kernel functions
    curl(l_i(x) g(z)):  a(curl(l_i g), curl(l_i g)) = |l_i|^2 |g'|^2 + m |l_i'|^2 |g|^2
    = (4h/15) [ |g'|^2 + kappa^2 |g|^2 ],  kappa^2 = 17.5 m / dx^2,  kappa^2 dz^2 = 17.5 eta.
    Reported: C_s = 1/lambda_min(P A) and lambda_max(P A).

(2) x-direction cross-check: CG_2 periodic, N_x cells, S_rho = M + rho h^2 K,
    (a) nodal (Jacobi) splitting alone, (b) vertex-star blocks + coarse CG_2(2h) (two-level),
    (c) hierarchical: coarse CG_2(2h) + fine bubbles (column patch), to reproduce C(rho) of the study.
"""
import numpy as np

# 1D quadratic Lagrange element on [0,1]: nodes 0, 1/2, 1 ; basis l0=(1-x)(1-2x), b=4x(1-x), l1=x(2x-1)
Mloc = (1/30.0)*np.array([[4, 2, -1],[2, 16, 2],[-1, 2, 4]], float)          # mass, times h
Kloc = (1/3.0)*np.array([[7, -8, 1],[-8, 16, -8],[1, -8, 7]], float)         # stiffness, times 1/h

def assemble_1d(N, periodic):
    # dof numbering: vertex l -> 2l, bubble l+1/2 -> 2l+1 ; N cells
    ndof = 2*N if periodic else 2*N+1
    M = np.zeros((ndof, ndof)); K = np.zeros((ndof, ndof))
    for c in range(N):
        idx = [2*c, 2*c+1, (2*c+2) % ndof if periodic else 2*c+2]
        for a in range(3):
            for b in range(3):
                M[idx[a], idx[b]] += Mloc[a, b]
                K[idx[a], idx[b]] += Kloc[a, b]
    return M, K

def asm_spectrum(A, patches, coarse=None):
    P = np.zeros_like(A)
    for p in patches:
        p = np.array(p)
        Aii = A[np.ix_(p, p)]
        P[np.ix_(p, p)] += np.linalg.inv(Aii)
    if coarse is not None:
        Ac = coarse.T @ A @ coarse
        P += coarse @ np.linalg.solve(Ac, coarse.T)
    ev = np.linalg.eigvals(P @ A).real
    ev = np.sort(ev)
    return ev[0], ev[-1]

print("=== (1) z-direction: 1D vertex-star ASM, no coarse space, Dirichlet both ends ===")
print("C_s = 1/lambda_min ; eta = m dz^2/dx^2 ; kappa^2 dz^2 = 17.5 eta ; rho_1 = eta N_z^2/pi^2")
print(f"{'N_z':>4} {'eta':>8} {'rho_1':>9} {'kappa*dz':>9} {'C_s':>9} {'lam_max':>8}   {'(N_z/pi)^2':>10}")
for Nz in (10, 20, 40, 80):
    h = 1.0/Nz
    M, K = assemble_1d(Nz, periodic=False)
    # Dirichlet: remove vertex 0 and vertex Nz
    keep = [d for d in range(2*Nz+1) if d not in (0, 2*Nz)]
    pos = {d: i for i, d in enumerate(keep)}
    Mh = h*M[np.ix_(keep, keep)]; Kh = K[np.ix_(keep, keep)]/h
    # patches: interior vertices l=1..Nz-1 : {2l, 2l-1, 2l+1}
    patches = [[pos[2*l], pos[2*l-1], pos[2*l+1]] for l in range(1, Nz)]
    for eta in (0.0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0):
        kappa2 = 17.5*eta/h**2
        A = Kh + kappa2*Mh
        lmin, lmax = asm_spectrum(A, patches)
        print(f"{Nz:4d} {eta:8.0e} {eta*Nz**2/np.pi**2:9.3g} {np.sqrt(kappa2)*h:9.3g} {1/lmin:9.2f} {lmax:8.2f}   {(Nz/np.pi)**2:10.1f}")

print()
print("=== (2) x-direction cross-checks, CG_2 periodic, S_rho = M + rho h^2 K ===")
Nx = 32
h = 1.0/Nx
M, K = assemble_1d(Nx, periodic=True)
ndof = 2*Nx
Mh = h*M; Kh = K/h
# coarse space CG_2(2h): prolongation by nodal interpolation (fine nodes are coarse nodes or coarse midpoints)
# coarse cell C spans fine cells 2C, 2C+1 : fine dofs 4C (vertex), 4C+1 (bubble), 4C+2 (vertex=coarse midpoint), 4C+3 (bubble), 4C+4 (vertex)
# coarse dofs: coarse vertex 2C -> cdof 2C, coarse midpoint -> cdof 2C+1 ; coarse ndof = 2*(Nx//2) = Nx
Ncx = Nx//2
Pc = np.zeros((ndof, Nx))
def q2_vals(xi):  # values of coarse basis (l0, b, l1) at local coordinate xi in [0,1]
    return np.array([(1-xi)*(1-2*xi), 4*xi*(1-xi), xi*(2*xi-1)])
for C in range(Ncx):
    cd = [2*C, 2*C+1, (2*C+2) % Nx]
    # fine nodes at xi = 0, 1/4, 1/2, 3/4, 1 of the coarse cell
    fine = [4*C, 4*C+1, 4*C+2, 4*C+3, (4*C+4) % ndof]
    for f, xi in zip(fine, (0, .25, .5, .75, 1.0)):
        if xi in (0, 1.0) and f != 4*C:   # shared vertex at xi=1 belongs to next coarse cell; set once
            continue
        Pc[f, :] = 0
        for a in range(3):
            Pc[f, cd[a]] += q2_vals(xi)[a]
# fix the vertices: Pc rows of coarse vertices (fine dof 4C) = identity on cdof 2C
for C in range(Ncx):
    Pc[4*C, :] = 0; Pc[4*C, 2*C] = 1.0
assert np.linalg.matrix_rank(Pc) == Nx

nodal = [[d] for d in range(ndof)]
stars = [[2*i, (2*i-1) % ndof, 2*i+1] for i in range(Nx)]          # vertex i with its two bubbles
bubbles = [[2*i+1] for i in range(Nx)]                               # column patch kernel part

print(f"{'rho':>8} | {'nodal alone: C':>14} | {'stars+coarse: lam_min':>22} {'lam_max':>8} | {'bubbles+coarse (hier.) C':>24} {'lam_max':>8}")
for rho in (0.0, 1e-3, 1e-2, 0.12, 1.0, 10.0, 1e4):
    S = Mh + rho*h*h*Kh
    lmin_n, _ = asm_spectrum(S, nodal)
    lmin_s, lmax_s = asm_spectrum(S, stars, coarse=Pc)
    lmin_b, lmax_b = asm_spectrum(S, bubbles, coarse=Pc)
    print(f"{rho:8.3g} | {1/lmin_n:14.2f} | {lmin_s:22.3f} {lmax_s:8.2f} | {1/lmin_b:24.2f} {lmax_b:8.2f}")

print()
print("=== (3) Q2 model: 1D Laplacian symbol, Jacobi vs overlapping 2-node (Vanka column) blocks ===")
th = np.linspace(np.pi/2, np.pi, 7)
jac = 1 - np.cos(th)                       # D^{-1}A symbol, A symbol 2-2cos
vanka = (4/3)*(1-np.cos(th))*(2+np.cos(th))
print("theta/pi   Jacobi PA   overlapping-pair PA")
for t, a, b in zip(th/np.pi, jac, vanka):
    print(f"{t:8.3f} {a:10.3f} {b:14.3f}")
wj = 2/(jac.min()+jac.max()); wv = 2/(vanka.min()+vanka.max())
print(f"optimal damping, smoothing factor on rough half: Jacobi {max(abs(1-wj*jac)):.3f}  pair-blocks {max(abs(1-wv*vanka)):.3f}")
