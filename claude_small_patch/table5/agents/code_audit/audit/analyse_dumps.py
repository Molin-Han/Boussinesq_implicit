"""Compare the operators dumped by lb_audit_copy.py (FS and MONO, nx=nz=4, serial). numpy only."""
import numpy as np, sys
fs = np.load(sys.argv[1]); mo = np.load(sys.argv[2])
rel = lambda A, B: np.linalg.norm(A - B) / np.linalg.norm(B)
n = fs['ises']; off = fs['ises_first']
print('field sizes (u_xz, u_y, b, p):', n, ' first global index per field:', off, ' use_rotation:', bool(fs['use_rotation']))
iv = np.arange(0, n[:3].sum()); ip = np.arange(n[:3].sum(), n.sum())
J, Jp = fs['J'], fs['Jp']
print('J == J(mono run):', np.array_equal(J, mo['J']), ' Jp == Jp(mono run):', np.array_equal(Jp, mo['Jp']))
D = Jp - J
print('Jp - J: nonzero only in (p,p) block:', np.abs(D[np.ix_(iv, iv)]).max() == 0 and np.abs(D[np.ix_(iv, ip)]).max() == 0 and np.abs(D[np.ix_(ip, iv)]).max() == 0,
      ' |J_pp|max =', np.abs(J[np.ix_(ip, ip)]).max())
A00, A01, A10, A11 = Jp[np.ix_(ip, ip)], Jp[np.ix_(ip, iv)], Jp[np.ix_(iv, ip)], Jp[np.ix_(iv, iv)]
blk = A00.copy()
for c in range(0, len(ip), 4): blk[c:c+4, c:c+4] = 0
print('delta*M_p block diagonal with 4x4 cell blocks:', np.abs(blk).max() == 0, ' symmetric:', np.allclose(A00, A00.T))
S = A11 - A10 @ np.linalg.solve(A00, A01)
print('Q1  FS aux vs exact Schur A11 - A10 A00^-1 A01 of Jp: rel Frobenius err = %.2e, max|diff|/max|S| = %.2e' % (rel(fs['aux'], S), np.abs(fs['aux'] - S).max() / np.abs(S).max()))
print('Q1  MONO aux vs Jp: rel err = %.2e' % rel(mo['aux'], Jp))
print('    fieldsplit_1 Pmat (pmat[1]) == A11 block of Jp:', np.array_equal(fs['fs1_P'], A11), '; fieldsplit_1 Amat type:', fs['fs1_Atype'])
bc = fs['bc_nodes']
print('    BC dofs (u_xz top/bottom):', len(bc), '; aux rows/cols identity on BC dofs:',
      np.array_equal(fs['aux'][bc][:, bc], np.eye(len(bc))) and np.abs(np.delete(fs['aux'][bc], bc, axis=1)).max() == 0 and np.abs(np.delete(fs['aux'][:, bc], bc, axis=0)).max() == 0)
x, y = fs['fs0_x'], fs['fs0_y']; z = np.linalg.solve(A00, x)
print('Q2  fieldsplit_0 vs (delta M_p)^-1 x: rel err %.2e ; vs mean-free (delta M_p)^-1 x: rel err %.2e ; mean(y) = %.1e' % (rel(y, z), rel(y, z - z.mean()), y.mean()))
print('    |A10 1| / |A10| =', np.abs(A10 @ np.ones(len(ip))).max() / np.abs(A10).max())
print('nullspaces: J', bool(fs['J_null']), ' Jp', bool(fs['Jp_null']), ' fieldsplit_0 Amat', bool(fs['fs0_A_null']), ' FS aux', bool(fs['aux_null']), ' MONO aux', bool(mo['aux_null']))
for tag, d in (('FS', fs), ('MONO', mo)):
    for l in range(3):
        print('  %-4s level %d: Amat type %-16s Amat is Pmat %-5s Amat is outer J %-5s Amat nullspace %-5s Pmat nullspace %-5s' % (
            tag, l, d[f'mg_L{l}_Atype'], bool(d[f'mg_L{l}_A_is_P']), bool(d[f'mg_L{l}_A_is_J']), bool(d[f'mg_L{l}_A_null']), bool(d[f'mg_L{l}_P_null'])))
print('Q3  level operators: FS level-l operator vs Schur complement (pressure eliminated) of MONO level-l operator')
for l in range(3):
    Pf, Pm = fs[f'mg_L{l}_P'], mo[f'mg_L{l}_P']; k = Pf.shape[0]
    Sm = Pm[:k, :k] - Pm[:k, k:] @ np.linalg.solve(Pm[k:, k:], Pm[k:, :k])
    print('    level %d: sizes FS %d, MONO %d (pressure %d): rel err %.2e' % (l, k, Pm.shape[0], Pm.shape[0] - k, rel(Pf, Sm)))
print('Q3  prolongations: MONO I_l restricted to (u,b) == FS I_l; MONO field coupling (ub<-p, p<-ub) zero')
for l in (1, 2):
    If, Im = fs[f'mg_L{l}_I'], mo[f'mg_L{l}_I']; kf, kc = If.shape
    print('    level %d<-%d: shape FS %s MONO %s: max|diff| %.1e (max|I| %.2f); |I_ub,p| %.1e |I_p,ub| %.1e; pressure block: max row sum %.3f min %.3f' % (
        l, l-1, If.shape, Im.shape, np.abs(Im[:kf, :kc] - If).max(), np.abs(If).max(), np.abs(Im[:kf, kc:]).max(), np.abs(Im[kf:, :kc]).max(),
        Im[kf:, kc:].sum(1).max(), Im[kf:, kc:].sum(1).min()))
print('Q4  patches: per level, MONO patch = FS patch + pressure dofs of the two columns; Schur(MONO patch) == FS patch')
for l in (1, 2):
    Kf, Km = fs[f'mg_L{l}_patches'], mo[f'mg_L{l}_patches']; k = Kf.shape[1]
    errs = []; lead = []
    for a, b in zip(Kf, Km):
        lead.append(np.abs(b[:k, :k] - a).max() / np.abs(a).max())
        Sb = b[:k, :k] - b[:k, k:] @ np.linalg.solve(b[k:, k:], b[k:, :k]); errs.append(rel(a, Sb))
    print('    level %d: %d patches, FS size %s, MONO size %s; max rel err Schur(MONO)-FS %.2e' % (
        l, len(Kf), sorted(set(fs[f'mg_L{l}_patchsize'])), sorted(set(mo[f'mg_L{l}_patchsize'])), max(errs)))
    for tag, d in (('FS', fs), ('MONO', mo)):
        print('      %-4s patch LU (preonly+lu, rcm) vs dense solve: max rel err %.1e ; cond(patch) min %.1e max %.1e' % (
            tag, d[f'mg_L{l}_patch_lu_relerr'].max(), d[f'mg_L{l}_patch_cond'].min(), d[f'mg_L{l}_patch_cond'].max()))
