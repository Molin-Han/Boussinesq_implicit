'''
Independent check (written without using the scripts in ../theory_kernel) of two statements
used in the write-up, for quadratic elements in x on a periodic mesh:

 (1) CG_2(h) = CG_2(2h) (+) span{cell bubbles of the fine mesh}
     (the fine x-bubbles span the complement of the coarse space), whereas for linear
     elements the complement of CG_1(2h) in CG_1(h) is spanned by the hats of the odd
     vertices, which are NOT bubbles of a single cell.

 (2) the constant of the splitting U = I_H U + (U - I_H U) (I_H = interpolation at the fine
     vertices by coarse quadratics) in the norm  S_rho(V) = |V|^2 + rho h^2 |V'|^2 :
         C(rho) = sup_U [ S_rho(I_H U) + S_rho(U - I_H U) ] / S_rho(U).
     rho -> 0   : the energy of the divergence-free fields is L2-like in x  (thin domain)
     rho -> inf : it is H1-like in x                                         (thick domain / large dt)

Hierarchical quadratic basis on [0, h]: vertex hats (1 - s), s and the bubble 4 s (1 - s).
'''
import numpy as np
from scipy.linalg import eigh


def element_matrices(h):
    # basis: phi0 = 1 - s, phi1 = s, phi2 = 4 s (1 - s),  s = x / h
    M = h * np.array([[1/3, 1/6, 1/3], [1/6, 1/3, 1/3], [1/3, 1/3, 8/15]])
    K = (1 / h) * np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 16/3]])
    return M, K


def assemble(N, h):
    '''periodic CG_2 on N cells: unknowns = N vertex values, then N bubble coefficients'''
    n = 2 * N
    M, K = np.zeros((n, n)), np.zeros((n, n))
    Me, Ke = element_matrices(h)
    for e in range(N):
        idx = [e, (e + 1) % N, N + e]
        M[np.ix_(idx, idx)] += Me
        K[np.ix_(idx, idx)] += Ke
    return M, K


def prolongation(Nc):
    '''CG_2(2h) -> CG_2(h) in the hierarchical bases above (coarse cell J = fine cells 2J, 2J+1)'''
    Nf = 2 * Nc
    P = np.zeros((2 * Nf, 2 * Nc))
    for J in range(Nc):
        vl, vr, b = J, (J + 1) % Nc, Nc + J
        # value of the coarse functions at the fine vertices 2J (s=0), 2J+1 (s=1/2)
        P[2 * J, vl] += 1.0
        P[2 * J + 1, vl] += 0.5
        P[2 * J + 1, vr] += 0.5
        P[2 * J + 1, b] += 1.0            # 4 s (1 - s) at s = 1/2
        # fine bubbles: coarse bubble restricted to a fine cell = its vertex interpolant + (1/4) fine bubble
        P[Nf + 2 * J, b] += 0.25
        P[Nf + 2 * J + 1, b] += 0.25
    return P


N = 16
h = 1.0
M, K = assemble(N, h)
P = prolongation(N // 2)
Mc, Kc = assemble(N // 2, 2 * h)
print('prolongation reproduces the coarse matrices: mass %.1e, stiffness %.1e'
      % (abs(P.T @ M @ P - Mc).max(), abs(P.T @ K @ P - Kc).max()))

# (1) dimension count: coarse space + fine bubbles
B = np.zeros((2 * N, N)); B[N:, :] = np.eye(N)                  # fine bubbles
rank2 = np.linalg.matrix_rank(np.hstack([P, B]))
print(f'k = 2: dim CG_2(h) = {2 * N}, rank [coarse | fine bubbles] = {rank2}  ->  spans the fine space: {rank2 == 2 * N}')
# linear elements: coarse hats + (no bubbles exist): complement = odd-vertex hats, each supported in TWO cells
print(f'k = 1: dim CG_1(h) = {N}, dim CG_1(2h) = {N // 2}; the {N // 2} missing functions are vertex hats (support = 2 cells): no single-cell function exists')

# (2) the splitting constant
R = np.zeros((N, 2 * N))                                          # coarse coefficients of I_H U
for J in range(N // 2):
    R[J, 2 * J] = 1.0                                             # coarse vertex J = fine vertex 2J
    R[N // 2 + J, 2 * J + 1] = 1.0                                # coarse bubble coefficient ...
    R[N // 2 + J, 2 * J] = -0.5                                   # ... = U(mid) - mean of the two ends
    R[N // 2 + J, (2 * J + 2) % N] = -0.5
IH = P @ R
Id = np.eye(2 * N)
# U - I_H U must be a combination of fine bubbles only
print('U - I_H U has no vertex component: %.1e' % abs((Id - IH)[:N, :]).max())
print('I_H is a projection onto the coarse space: %.1e' % abs(IH @ IH - IH).max())
print()
print('   rho        C(rho)      1/C(rho)')
for rho in [0.0, 1e-3, 1e-2, 0.12, 1.0, 2.0, 10.0, 1e3, np.inf]:
    if np.isinf(rho):
        S = K + 1e-12 * M
    else:
        S = M + rho * h * h * K
    T = IH.T @ S @ IH + (Id - IH).T @ S @ (Id - IH)
    C = eigh(T, S, eigvals_only=True)[-1]
    print(f'{rho:8.3g}   {C:9.4f}   {1 / C:9.4f}')

# several levels: hierarchical splitting U = U_0 + sum_l (bubbles of level l), L2-like energy (rho = 0) and H1-like
print()
print('number of levels, splitting into coarsest space + bubbles of every finer level (N = 64 on the finest level):')
for rho in [0.0, 1.0, np.inf]:
    out = []
    for nlev in [2, 3, 4, 5]:
        Nf = 64
        Mf, Kf = assemble(Nf, 1.0)
        S = (Kf + 1e-12 * Mf) if np.isinf(rho) else (Mf + rho * Kf)
        # pieces: Q_l = I_l - I_{l-1} expressed on the finest level, I_finest = identity
        Pl = np.eye(2 * Nf)          # prolongation level l -> finest
        Il_prev = np.eye(2 * Nf)     # I_l on the finest level (starts with the identity)
        T = np.zeros_like(S)
        n = Nf
        for lev in range(nlev - 1):
            Pc = prolongation(n // 2)
            Rc = np.zeros((n, 2 * n))
            for J in range(n // 2):
                Rc[J, 2 * J] = 1.0
                Rc[n // 2 + J, 2 * J + 1] = 1.0
                Rc[n // 2 + J, 2 * J] = -0.5
                Rc[n // 2 + J, (2 * J + 2) % n] = -0.5
            # interpolation of a level-(current) function to the next coarser level, lifted to the finest level
            # current-level coefficients of I_cur U:  Ccur (n_cur x n_finest)
            if lev == 0:
                Ccur = np.eye(2 * Nf)
            Ccoarse = Rc @ Ccur                       # coarse coefficients
            Inext = (Pl @ Pc) @ Ccoarse               # I_{coarser} on the finest level
            Q = Il_prev - Inext                       # bubbles of the current level
            T += Q.T @ S @ Q
            Il_prev, Pl, Ccur, n = Inext, Pl @ Pc, Ccoarse, n // 2
        T += Il_prev.T @ S @ Il_prev                  # coarsest component
        out.append(eigh(T, S, eigvals_only=True)[-1])
    print(f'   rho = {rho:6.3g}:  C for 2, 3, 4, 5 levels = ' + ', '.join(f'{c:8.3f}' for c in out))
