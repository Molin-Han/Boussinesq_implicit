"""numerics agent: why do FS and MONO-P give identical velocities with the GMRES(6) level smoother?
Runs the set-up of fs_vs_mono.py for two solvers only (FGMRES outer, GMRES(6) on MG levels) and compares
the finest-level smoothers in isolation (PETSc objects) and a numpy GMRES on the same operators."""
import sys, os
sys.argv = [sys.argv[0], '--nx', os.environ.get('NX', '4'), '--nz', os.environ.get('NZ', '4'),
            '--counts_only', '--variants', 'gmres', '--methods', 'FS,MONO-P', '--out', 'gmres_diag_counts.json']
sys.dont_write_bytecode = True
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fs_vs_mono.py')).read()
_src = _src.replace("if args.counts_only:\n    with open(OUT, 'w') as f:\n        json.dump(results, f, indent=1)\n    sys.exit(0)", "")
_src = _src.split("# Operators actually used on the MG levels")[0]       # stop after part (c)
exec(compile(_src, 'fs_vs_mono.py', 'exec'))

sF = kept[('gmres', 'FS')]
Jmat, Jpmat = sF.snes.ksp.getOperators()
J_d = to_dense(Jmat); Jp_d = to_dense(Jpmat)
sM = kept[('gmres', 'MONO-P')]
mgF, ctxF = inner_mg(sF, 'FS')
mgM, ctxM = inner_mg(sM, 'MONO-P')
kF = mgF.getMGSmoother(mgF.getMGLevels() - 1)
kM = mgM.getMGSmoother(mgM.getMGLevels() - 1)
pcF = kF.getPC(); pcM = kM.getPC()
A00 = Jp_d[np.ix_(IP, IP)]; A01 = Jp_d[np.ix_(IP, IV)]; A10 = Jp_d[np.ix_(IV, IP)]; A11 = Jp_d[np.ix_(IV, IV)]
S_d = A11 - A10 @ np.linalg.solve(A00, A01)
rng = np.random.default_rng(7)
r = rng.standard_normal(NTOT); r[BCDOFS] = 0.0
rt = r[IV] - A10 @ np.linalg.solve(A00, r[IP])

def pcapply(pc, Amat, v, solver):
    xv, yv = Amat.createVecs()
    xv.setArray(v.copy()); yv.set(0.0)
    with solver_context(solver):
        pc.apply(xv, yv)
    return yv.getArray().copy()

AF = kF.getOperators()[1]; AM = kM.getOperators()[1]
MinvF = lambda v: pcapply(pcF, AF, v, sF)
MinvM = lambda v: pcapply(pcM, AM, v, sM)

# 1) the PETSc smoother objects, zero initial guess
def ksp_apply(k, Amat, b, solver):
    xv, bv = Amat.createVecs()
    bv.setArray(b.copy()); xv.set(0.0)
    with solver_context(solver):
        k.solve(bv, xv)
    return xv.getArray().copy()
xF = ksp_apply(kF, AF, rt, sF)
xM = ksp_apply(kM, AM, r, sM)
print("PETSc finest GMRES(6) smoother alone: ||xM_v - xF||/||xF|| =", np.linalg.norm(xM[:len(IV)] - xF) / np.linalg.norm(xF))
print("   error of xF vs S^{-1} rt:", np.linalg.norm(xF - np.linalg.solve(S_d, rt)) / np.linalg.norm(np.linalg.solve(S_d, rt)))

# 2) numpy left-preconditioned GMRES(k), MGS, zero guess
def gmres(Aop, Minv, b, k, proj=None):
    z0 = Minv(b)
    if proj: z0 = proj(z0)
    beta = np.linalg.norm(z0)
    V = [z0 / beta]; H = np.zeros((k + 1, k))
    for j in range(k):
        w = Minv(Aop(V[j]))
        if proj: w = proj(w)
        for i in range(j + 1):
            H[i, j] = np.dot(w, V[i]); w = w - H[i, j] * V[i]
        H[j + 1, j] = np.linalg.norm(w); V.append(w / H[j + 1, j])
    e1 = np.zeros(k + 1); e1[0] = beta
    y = np.linalg.lstsq(H, e1, rcond=None)[0]
    res = [np.linalg.norm(e1[:i+2] - H[:i+2, :i+1] @ np.linalg.lstsq(H[:i+2, :i+1], e1[:i+2], rcond=None)[0]) / beta for i in range(k)]
    return sum(yi * vi for yi, vi in zip(y, V[:k])), res
for k in [1, 2, 3, 6]:
    gF, resF = gmres(lambda v: S_d @ v, MinvF, rt, k)
    gM, resM = gmres(lambda v: Jp_d @ v, MinvM, r, k)
    print(f"numpy GMRES({k}): ||gM_v - gF||/||gF|| = {np.linalg.norm(gM[:len(IV)] - gF) / np.linalg.norm(gF):.3e}; "
          f"rel prec. residual FS {resF[-1]:.3e}, MONO-P {resM[-1]:.3e}")
# 3) Krylov space structure: the monolithic preconditioned residual z=M^{-1}r; is its pressure part a function of its velocity part?
z = MinvM(r)
print("z_p - [2 A00^{-1} r_p - A00^{-1}A01 z_v] rel:", np.linalg.norm(z[len(IV):] - (2*np.linalg.solve(A00, r[IP]) - np.linalg.solve(A00, A01 @ z[:len(IV)]))) / np.linalg.norm(z[len(IV):]))
print("norms: ||z_v||", np.linalg.norm(z[:len(IV)]), " ||z_p||", np.linalg.norm(z[len(IV):]))
