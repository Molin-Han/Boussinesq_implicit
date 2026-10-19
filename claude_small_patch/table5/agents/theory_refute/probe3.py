# Dense checks (small mesh): C3/C4 algebra on the finest level, spectra, and what the actual MONO cycle approximates.
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'probe.py')).read().split('solvers = {}')[0]
exec(src)
import numpy as np
nrm = np.linalg.norm
T = args.tag
ONLY = os.environ.get('P3_ONLY', '')

sv = {}
for name in ('fs', 'mono'):
    s = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=CONF[name], appctx=appctx, options_prefix=f'{name}_')
    Un.assign(IC); Unp1.assign(Un); s.solve()
    print(f'[{T}] {name} first-step its {s.snes.ksp.getIterationNumber()}')
    sv[name] = s
sF, sM = sv['fs'], sv['mono']
mgF = inner_mg(sF, 'fs'); mgM = inner_mg(sM, 'mono')
lev = mgM.getMGLevels() - 1
kM = mgM.getMGSmoother(lev); kF = mgF.getMGSmoother(lev)
dense = lambda Mt: Mt.convert('dense').getDenseArray().copy()
Ad = dense(kM.getOperators()[1])          # shifted A_delta (MONO Pmat)
Jd = dense(sM.snes.ksp.getOperators()[0]) # unshifted J (outer Amat == MONO finest Amat)
Sd = dense(kF.getOperators()[1])          # S_aux (FS Pmat)
nv = Sd.shape[0]; N = Ad.shape[0]; npr = N - nv
G = Ad[:nv, nv:]; D = Ad[nv:, :nv]; dM = Ad[nv:, nv:]
print(f'[{T}] N={N} nv={nv} np={npr}; |J - A_delta| only in pp block: {nrm((Jd-Ad)[:nv,:])+nrm((Jd-Ad)[:, :nv]):.2e}; |Jpp|={nrm(Jd[nv:,nv:]):.2e}')
Schur = Ad[:nv, :nv] - G @ np.linalg.solve(dM, D)
print(f'[{T}] |S_aux - Schur(A_delta)|/|S_aux| = {nrm(Sd-Schur)/nrm(Sd):.2e}')
from firedrake.dmhooks import get_function_space
def patches(ksp):
    pc_ = ksp.getPC(); ctx_ = pc_.getPythonContext(); V_ = get_function_space(pc_.getDM())
    with solver_ctx(sM):
        ises_ = ctx_.get_patches(V_)
    lg = V_.dof_dset.lgmap
    return [lg.applyIS(i).getIndices().copy() for i in ises_]
isM = patches(kM); isF = patches(kF)
print(f'[{T}] patches: MONO {len(isM)} sizes {sorted(set(len(i) for i in isM))}; FS {len(isF)} sizes {sorted(set(len(i) for i in isF))}; velocity part of MONO patch == FS patch: {sum(set(a[a < nv].tolist()) == set(b_.tolist()) for a, b_ in zip(isM, isF))}/{len(isF)}')
mult = np.zeros(N)
for idx in isM: mult[idx] += 1
print(f'[{T}] multiplicity of dofs in MONO patches: velocity {sorted(set(mult[:nv].astype(int).tolist()))}, pressure {sorted(set(mult[nv:].astype(int).tolist()))}')
def asm_dense(Amat, ises):
    B = np.zeros_like(Amat)
    for idx in ises:
        B[np.ix_(idx, idx)] += np.linalg.inv(Amat[np.ix_(idx, idx)])
    return B
BM = asm_dense(Ad, isM); BS = asm_dense(Sd, isF)
TM = BM @ Ad; TS = BS @ Sd; TJ = BM @ Jd
print(f'[{T}] C4: |TM[v,p]|/|TM| = {nrm(TM[:nv, nv:])/nrm(TM):.2e}; |TM[v,v]-TS|/|TS| = {nrm(TM[:nv,:nv]-TS)/nrm(TS):.2e}; |TM[p,p]-2I|/|2I| = {nrm(TM[nv:,nv:]-2*np.eye(npr))/nrm(2*np.eye(npr)):.2e}')
L = np.block([[np.eye(nv), np.zeros((nv, npr))], [np.linalg.solve(dM, D), np.eye(npr)]])   # (e_v, d=e_p+(dM)^-1 D e_v)
TMd = L @ TM @ np.linalg.inv(L)
print(f'[{T}] (e_v,d) coords: |offdiag|/|T| = {(nrm(TMd[:nv,nv:])+nrm(TMd[nv:,:nv]))/nrm(TMd):.2e}; |X| = |TM[p,v]| = {nrm(TM[nv:,:nv]):.2e} vs |TM[v,v]| = {nrm(TM[:nv,:nv]):.2e}')
for nm, Tm in ((() if ONLY else (('ASM_S S', TS), ('ASM_M A_delta', TM), ('ASM_M J (actual MONO finest)', TJ)))):
    ev = np.linalg.eigvals(Tm)
    n2 = np.sum(np.abs(ev - 2.0) < 1e-6)
    H = 0.5*(Tm + Tm.T); fov = np.linalg.eigvalsh(H)
    print(f'[{T}] spec {nm}: Re in [{ev.real.min():.4f}, {ev.real.max():.4f}], max|Im|={np.abs(ev.imag).max():.3e}, #(|lam-2|<1e-6)={n2}/{len(ev)}, #(|lam|<1e-8)={np.sum(np.abs(ev)<1e-8)}; Euclidean field of values Re in [{fov.min():.3e}, {fov.max():.3e}]')
# What does each preconditioner approximate?  first-step RHS b (= -F(U^n)), exact shifted and unshifted solutions
b = sM.snes.ksp.getRhs().getArray().copy() if sM.snes.ksp.getRhs() is not None else None
Un.assign(IC); Unp1.assign(Un)
Fc = assemble(eqn, bcs=None)
with Fc.dat.vec_ro as fv:
    b = -fv.getArray().copy()
bcn = np.concatenate([W._ises[0].indices[bc.nodes] for bc in bcs]); b[bcn] = 0.0
iv = np.arange(nv); ip = np.arange(nv, N)
print(f'[{T}] RHS b: |b_v|={nrm(b[iv]):.3e} |b_p|={nrm(b[ip]):.3e}')
xd = np.linalg.solve(Ad, b)
# unshifted: bordered system with zero-mean pressure (J has constant-pressure null space)
nvec = np.zeros(N); nvec[ip] = 1.0/np.sqrt(npr)
Jb = np.block([[Jd, nvec[:, None]], [nvec[None, :], np.zeros((1, 1))]])
x0 = np.linalg.solve(Jb, np.concatenate([b, [0.0]]))[:N]
def pm0(xx):
    y_ = xx.copy(); y_[ip] -= y_[ip].mean(); return y_
def applyP(s, r):
    xv = sM.snes.ksp.getOperators()[0].createVecLeft(); xv.setArray(r); yv = xv.duplicate()
    with solver_ctx(s):
        s.snes.ksp.getPC().apply(xv, yv)
    return yv.getArray().copy()
xF = applyP(sF, b); xM = applyP(sM, b)
def err(xx, ref):
    xx = pm0(xx); ref = pm0(ref)
    return nrm(xx[iv]-ref[iv])/nrm(ref[iv]), nrm(xx[ip]-ref[ip])/nrm(ref[ip])
for nm, xx in (('A_delta^-1 b (exact shifted)', xd), ('P_FS^-1 b', xF), ('P_MONO^-1 b', xM)):
    ev_, ep_ = err(xx, x0); evd, epd = err(xx, xd)
    print(f'[{T}] {nm:30s}: vs unshifted x0: vel {ev_:.3e} pres {ep_:.3e} | vs shifted x_delta: vel {evd:.3e} pres {epd:.3e} | |J x - b|/|b| = {nrm(Jd@xx-b)/nrm(b):.3e}')
# preconditioned-operator spectra (dense) for outer GMRES: P^-1 J restricted to the complement of the null space
def dense_pc(s):
    Pm = np.zeros((N, N))
    for j in range(N):
        e = np.zeros(N); e[j] = 1.0; Pm[:, j] = applyP(s, e)
    return Pm
if N <= 1500:
    for nm, s in [(a_, b_) for a_, b_ in (('FS', sF), ('MONO', sM)) if (not ONLY or a_ == ONLY)]:
        Pinv = dense_pc(s)
        K = Pinv @ Jd
        ev = np.linalg.eigvals(K)
        ev = ev[np.abs(ev) > 1e-10]
        d1 = np.abs(ev - 1.0)
        print(f'[{T}] eig(P^-1 J) {nm}: min|lam|={np.abs(ev).min():.3e}, max|lam|={np.abs(ev).max():.3e}, #|lam-1|>0.1: {np.sum(d1>0.1)}, #|lam-1|>0.01: {np.sum(d1>0.01)}, #|lam-1|>1e-3: {np.sum(d1>1e-3)} of {len(ev)}')
    Pex = np.linalg.inv(Ad); K = Pex @ Jd; ev = np.linalg.eigvals(K); ev = ev[np.abs(ev) > 1e-10]; d1 = np.abs(ev-1.0)
    print(f'[{T}] eig(A_delta^-1 J) exact: min|lam|={np.abs(ev).min():.3e}, #|lam-1|>0.1: {np.sum(d1>0.1)}, #|lam-1|>0.01: {np.sum(d1>0.01)}, #|lam-1|>1e-3: {np.sum(d1>1e-3)} of {len(ev)}')
