# Level-wise checks of claims C3/C4: finest-level ASM of MONO vs FS, Schur identity, exact solves.
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'probe.py')).read().split('solvers = {}')[0]
exec(src)
import numpy as np
nrm = np.linalg.norm

sv = {}
for name in ('fs', 'mono'):
    s = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=CONF[name], appctx=appctx,
                                   options_prefix=f'{name}_')
    Un.assign(IC); Unp1.assign(Un); s.solve()
    print(f'[{args.tag}] {name} first-step its {s.snes.ksp.getIterationNumber()}')
    sv[name] = s
sF, sM = sv['fs'], sv['mono']
mgF = inner_mg(sF, 'fs'); mgM = inner_mg(sM, 'mono')
L = mgM.getMGLevels()
iv = np.concatenate([W._ises[i].indices for i in range(3)]); ip = W._ises[3].indices
rng = np.random.RandomState(7)

def mat_dense(Mt):
    return Mt.convert('dense').getDenseArray().copy()

for lev in [L-1]:
    kM = mgM.getMGSmoother(lev); kF = mgF.getMGSmoother(lev)
    AM = kM.getOperators()[1]; SF = kF.getOperators()[1]
    print(f'[{args.tag}]   MONO smoother Amat is outer J: {kM.getOperators()[0].handle == sM.snes.ksp.getOperators()[0].handle}; FS smoother Amat type: {kF.getOperators()[0].getType()}')
    print(f'[{args.tag}] level {lev}: MONO op size {AM.getSize()}, FS op size {SF.getSize()}, types {AM.getType()} {SF.getType()}')
    print(f'[{args.tag}]   MONO op nullspace attached: {AM.getNullSpace().handle != 0}; FS op nullspace: {SF.getNullSpace().handle != 0}')
    Ad = mat_dense(AM); Sd = mat_dense(SF)
    nv = Sd.shape[0]
    # In serial the MONO ordering is [u_xz, u_y, b, p] and the first nv entries are the FS ordering
    Avv = Ad[:nv, :nv]; G = Ad[:nv, nv:]; Dm = Ad[nv:, :nv]; dM = Ad[nv:, nv:]
    Schur = Avv - G @ np.linalg.solve(dM, Dm)
    print(f'[{args.tag}]   |S_aux - Schur(A_delta)|_F/|S_aux|_F = {nrm(Sd - Schur)/nrm(Sd):.3e}')
    print(f'[{args.tag}]   |Avv|_F={nrm(Avv):.3e} |G|_F={nrm(G):.3e} |D|_F={nrm(Dm):.3e} |dM|_F={nrm(dM):.3e}  |S|_F={nrm(Sd):.3e}')
    # finest ASM PCs
    asmM = kM.getPC().getPythonContext().asmpc; asmF = kF.getPC().getPythonContext().asmpc
    isM = asmM.getASMLocalSubdomains()[0]; isF = asmF.getASMLocalSubdomains()[0]
    print(f'[{args.tag}]   #patches MONO {len(isM)} FS {len(isF)}; sizes MONO {sorted(set(i.getSize() for i in isM))} FS {sorted(set(i.getSize() for i in isF))}')
    # compare patch index sets: MONO patch minus pressure dofs == FS patch?
    same = 0
    for a_, b_ in zip(isM, isF):
        sa = set(a_.getIndices().tolist()); sb = set(b_.getIndices().tolist())
        sav = set(i for i in sa if i < nv); sap = set(i for i in sa if i >= nv)
        same += (sav == sb)
    print(f'[{args.tag}]   patches whose velocity part equals the FS patch: {same}/{len(isF)}')
    # pressure dofs per patch and multiplicity
    mult = np.zeros(Ad.shape[0])
    for a_ in isM: mult[a_.getIndices()] += 1
    print(f'[{args.tag}]   MONO multiplicity: velocity dofs {sorted(set(mult[:nv].astype(int)))}, pressure dofs {sorted(set(mult[nv:].astype(int)))}')
    def apply_pc(pc, vec_np, Mt):
        xv = Mt.createVecRight(); xv.setArray(vec_np); yv = Mt.createVecLeft()
        pc.apply(xv, yv); return yv.getArray().copy()
    with solver_ctx(sM):
        pass
    for trial in range(2):
        rho = rng.randn(Ad.shape[0])
        rS = rho[:nv] - G @ np.linalg.solve(dM, rho[nv:])
        yM = apply_pc(asmM, rho, AM); yF = apply_pc(asmF, rS, SF)
        print(f'[{args.tag}]   ASM test {trial}: |ASM_M(rho)_v - ASM_S(r_S)|/|ASM_S(r_S)| = {nrm(yM[:nv]-yF)/nrm(yF):.3e}')
        # local patch check on patch 0: direct dense solves
        idxM = isM[0].getIndices(); idxF = isF[0].getIndices()
        cM = np.linalg.solve(Ad[np.ix_(idxM, idxM)], rho[idxM])
        cF = np.linalg.solve(Sd[np.ix_(idxF, idxF)], rS[idxF])
        posM = {g: k for k, g in enumerate(idxM)}
        cMv = np.array([cM[posM[g]] for g in idxF])
        print(f'[{args.tag}]   dense patch-0 solve: |cM_v - cF|/|cF| = {nrm(cMv-cF)/nrm(cF):.3e};  cond(A_ii)={np.linalg.cond(Ad[np.ix_(idxM, idxM)]):.2e} cond(S_ii)={np.linalg.cond(Sd[np.ix_(idxF, idxF)]):.2e}')
    # block-triangular structure of ASM_M A_delta (C4): velocity part depends only on e_v
    e = rng.randn(Ad.shape[0]); e2 = e.copy(); e2[nv:] = rng.randn(Ad.shape[0]-nv)
    t1 = apply_pc(asmM, Ad @ e, AM); t2 = apply_pc(asmM, Ad @ e2, AM)
    tS = apply_pc(asmF, Sd @ e[:nv], SF)
    print(f'[{args.tag}]   C4: |(T e)_v - (T e2)_v|/|(T e)_v| = {nrm(t1[:nv]-t2[:nv])/nrm(t1[:nv]):.3e}; |(T e)_v - ASM_S S e_v|/| | = {nrm(t1[:nv]-tS)/nrm(tS):.3e}')
    # pressure block 2I: e with e_v=0
    e3 = np.zeros(Ad.shape[0]); e3[nv:] = rng.randn(Ad.shape[0]-nv)
    t3 = apply_pc(asmM, Ad @ e3, AM)
    print(f'[{args.tag}]   C4: e_v=0 -> |(T e)_v|={nrm(t3[:nv]):.3e}, |(T e)_p - 2 e_p|/|2e_p| = {nrm(t3[nv:]-2*e3[nv:])/nrm(2*e3[nv:]):.3e}')
    # deviation d = e_p + (dM)^{-1} D e_v evolves as 2d?
    d0 = e2[nv:] + np.linalg.solve(dM, Dm @ e2[:nv])
    Te = apply_pc(asmM, Ad @ e2, AM)
    d1 = Te[nv:] + np.linalg.solve(dM, Dm @ Te[:nv])
    print(f'[{args.tag}]   (e_v,d) coords: |d(Te) - 2 d(e)|/|2d(e)| = {nrm(d1-2*d0)/nrm(2*d0):.3e}')
    np.save(os.path.join(HERE, f'dense_{args.tag}_A.npy'), Ad); np.save(os.path.join(HERE, f'dense_{args.tag}_S.npy'), Sd)
