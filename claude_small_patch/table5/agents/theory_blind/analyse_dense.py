"""Numerical check of the block analysis (theory_blind).  numpy only.
Loads runs/D_{FS,MONOnoamat,MONO}_n4/dense.npz (P^{-1}, A, A_delta) and runs/P_n4/pcols.npz."""
import numpy as np, scipy.linalg as sla, sys
H = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_blind/runs/'
tag = sys.argv[1] if len(sys.argv) > 1 else 'n4'
R = {k: np.load(H + f'D_{k}_{tag}/dense.npz') for k in ['FS', 'MONOnoamat', 'MONO']}
pc = np.load(H + f'P_{tag}/pcols.npz')
A, Ad = R['FS']['A'], R['FS']['Ad']
for k in R:
    assert np.allclose(R[k]['A'], A, rtol=0, atol=1e-12 * abs(A).max()) and np.allclose(R[k]['Ad'], Ad, rtol=0, atol=1e-12 * abs(Ad).max())
sizes = R['FS']['sizes']; delta = float(R['FS']['shift']); dt = float(R['FS']['dt'])
n = A.shape[0]; nv = int(sizes[:3].sum()); iv = np.arange(nv); ip = np.arange(nv, n)
print(f'n={n} sizes(u_xz,u_y,b,p)={list(sizes)} delta={delta} dt={dt}')
# blocks
D = A[np.ix_(ip, iv)]; Mp = (Ad - A)[np.ix_(ip, ip)] / delta; Avp = A[np.ix_(iv, ip)]
print('check A_vp = -dt D^T :', np.abs(Avp + dt * D.T).max() / np.abs(Avp).max())
print('check A_pp = 0       :', np.abs(A[np.ix_(ip, ip)]).max())
S = Ad[np.ix_(iv, iv)] - Avp @ np.linalg.solve(Ad[np.ix_(ip, ip)], D)
# T = [[I,0],[(1/delta) Mp^{-1} D, I]]
Tvp = np.linalg.solve(Mp, D) / delta
T = np.eye(n); T[np.ix_(ip, iv)] = Tvp
Ti = np.eye(n); Ti[np.ix_(ip, iv)] = -Tvp
# projector removing the constant-pressure direction (the null vector of A), used on both sides
c = np.zeros(n); c[ip] = 1.0; c /= np.linalg.norm(c); Pi = np.eye(n) - np.outer(c, c)
def blocks(M):
    M = Pi @ M @ Pi
    return {'vv': M[np.ix_(iv, iv)], 'vp': M[np.ix_(iv, ip)], 'pv': M[np.ix_(ip, iv)], 'pp': M[np.ix_(ip, ip)]}
nrm = lambda X: np.linalg.norm(X, 2)
TE = {}
for k in R:
    E = np.eye(n) - R[k]['Pinv'] @ Ad
    TE[k] = blocks(T @ E @ Ti)
    b = TE[k]
    print(f'[{k}] T(I-P^-1 A_delta)T^-1 block 2-norms: vv {nrm(b["vv"]):.3e}  v-eta {nrm(b["vp"]):.3e}  eta-v {nrm(b["pv"]):.3e}  eta-eta {nrm(b["pp"]):.3e}'
          f'   rho(vv) {max(abs(np.linalg.eigvals(b["vv"]))):.3e}')
print('||TE_vv(MONOnoamat) - TE_vv(FS)|| / ||TE_vv(FS)|| =', nrm(TE['MONOnoamat']['vv'] - TE['FS']['vv']) / nrm(TE['FS']['vv']))
# spectra of P^{-1} A (drop the single zero eigenvalue of the pressure null space)
def spec(M):
    ev = np.linalg.eigvals(M); ev = ev[np.argsort(abs(ev))]
    return ev[1:]
G = np.linalg.solve(Ad, A)
specs = {'ideal A_delta^-1 A': spec(G)}
for k in R: specs[k] = spec(R[k]['Pinv'] @ A)
for k, ev in specs.items():
    out = ev[abs(ev - 1) > 0.1]
    print(f'{k:20s}: #eig with |lam-1|>0.1: {len(out):4d};  min|lam| {abs(ev).min():.3e};  max|lam| {abs(ev).max():.3f};  max|Im| {abs(ev.imag).max():.2e}')
    print('      smallest 12 |lam|:', np.array2string(np.sort(abs(ev))[:12], precision=4, max_line_width=200))
# classify the pressure eigenvalues of the ideal operator: barotropic fraction of eigenvectors
w, V = np.linalg.eig(G)
xc = np.round(pc['xc'], 3); xp = np.round(pc['xp'], 3)
keys = {}
gid = np.array([keys.setdefault((a, b), len(keys)) for a, b in zip(xc, xp)])   # (column, x-node) groups
def baro_frac(p):
    m = np.zeros(len(keys), dtype=complex); np.add.at(m, gid, p); cnt = np.bincount(gid); m /= cnt
    return np.linalg.norm(m[gid])**2 / max(np.linalg.norm(p)**2, 1e-300)
sel = np.where((abs(w) > 1e-10) & (abs(w - 1) > 1e-6))[0]
sel = sel[np.argsort(w[sel].real)]
print('ideal A_delta^-1 A: eigenvalues != 1 (pressure modes), with barotropic fraction of the pressure part:')
rows = [(w[i].real, baro_frac(V[ip, i])) for i in sel]
for lam, bf in rows[:40]:
    print(f'   mu={lam:.4e}   baro={bf:.3f}')
print('   ... total', len(rows), 'pressure eigenvalues; mu range', rows[0][0], rows[-1][0])
