"""numerics agent: block-by-block (per field) versions of checks (a) and (b) of fs_vs_mono.py, nx=nz=4.
(a) ||S_aux - (A11 - A10 A00^{-1} A01)|| per field block (u_xz, u_y, b), also on the coarse MG levels;
(b) per-field relative differences ||x_f^M - x_f^FS|| / ||x_f^FS|| for f in (u_xz, u_y, b, p mod const)."""
import sys, os
sys.argv = [sys.argv[0], '--nx', '4', '--nz', '4', '--counts_only', '--variants', 'rich,chebfix,chebest',
            '--methods', 'FS,MONO,MONO-P', '--sides', 'left', '--out', 'blockwise_counts.json']
sys.dont_write_bytecode = True
_src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fs_vs_mono.py')).read()
_src = _src.split("# Operators actually used on the MG levels")[0]
_src = _src.replace("if args.counts_only:\n    with open(OUT, 'w') as f:\n        json.dump(results, f, indent=1)\n    sys.exit(0)", "")
exec(compile(_src, 'fs_vs_mono.py', 'exec'))

s0 = kept[('chebest', 'FS')]
Jmat, Jpmat = s0.snes.ksp.getOperators()
J_d = to_dense(Jmat); Jp_d = to_dense(Jpmat)
A00 = Jp_d[np.ix_(IP, IP)]; A01 = Jp_d[np.ix_(IP, IV)]; A10 = Jp_d[np.ix_(IV, IP)]; A11 = Jp_d[np.ix_(IV, IV)]
GD = A10 @ np.linalg.solve(A00, A01)
S_d = A11 - GD
mg_fs, ctx_fs = inner_mg(s0, 'FS')
Saux = to_dense(ctx_fs.P.petscmat)
names = ['u_xz', 'u_y', 'b']
off = np.cumsum([0] + [len(f) for f in fis[:3]])
print("(a) fine level, per block: ||Saux_ij - S_ij||_F / ||S_ij||_F   [||S_ij||_F]")
for i in range(3):
    for j in range(3):
        bi = slice(off[i], off[i+1]); bj = slice(off[j], off[j+1])
        nS = np.linalg.norm(S_d[bi, bj]); nd = np.linalg.norm(Saux[bi, bj] - S_d[bi, bj])
        print(f"   {names[i]:5s}-{names[j]:5s}: {nd / nS if nS > 0 else float('nan'):.3e}   [{nS:.3e}]  abs diff {nd:.2e}")
bx = slice(off[0], off[1])
print(f"   u_xz-u_xz: ||diff||/||A10 A00^-1 A01||_F = {np.linalg.norm(Saux[bx, bx] - S_d[bx, bx]) / np.linalg.norm(GD[bx, bx]):.3e}"
      f"  (||grad-div part||_F={np.linalg.norm(GD[bx, bx]):.3e}, ||A11_xx||_F={np.linalg.norm(A11[bx, bx]):.3e})")
# coarse levels: Schur complement of the MONO-P level operator vs FS level operator, u_xz block only
mgM, ctxM = inner_mg(kept[('chebest', 'MONO-P')], 'MONO-P')
for l in range(mg_fs.getMGLevels()):
    kf = mg_fs.getMGSmoother(l) if l > 0 else mg_fs.getMGCoarseSolve()
    km = mgM.getMGSmoother(l) if l > 0 else mgM.getMGCoarseSolve()
    Sf = to_dense(kf.getOperators()[1]); Am = to_dense(km.getOperators()[1])
    nv = Sf.shape[0]
    gd = Am[:nv, nv:] @ np.linalg.solve(Am[nv:, nv:], Am[nv:, :nv])
    Sl = Am[:nv, :nv] - gd
    nxz = (Am.shape[0] - nv) // 4 * 2 + 0   # placeholder
    print(f"   level {l}: ||S_aux,l - S_l||_F/||S_l||_F = {np.linalg.norm(Sf - Sl) / np.linalg.norm(Sl):.3e}; relative to ||grad-div part||_F: {np.linalg.norm(Sf - Sl) / np.linalg.norm(gd):.3e}")

# (b) per field
rng = np.random.default_rng(20261019)
r_rand = rng.standard_normal(NTOT); r_rand[BCDOFS] = 0.0
Fc = assemble(eqn)
with Fc.dat.vec_ro as fv:
    r_res = fv.getArray().copy()
r_res[BCDOFS] = 0.0
def apply_pc(solver, r):
    pc = solver.snes.ksp.getPC()
    xv, yv = solver.snes.ksp.getOperators()[0].createVecs()
    xv.setArray(r.copy()); yv.set(0.0)
    with solver_context(solver):
        pc.apply(xv, yv)
    return yv.getArray().copy()
fields = [('u_xz', fis[0]), ('u_y', fis[1]), ('b', fis[2]), ('p', fis[3])]
print("(b) per-field relative difference vs FS: ||x_f - x_f^FS|| / ||x_f^FS||  (p: constant removed)")
for rname, r in [('random', r_rand), ('residual', r_res)]:
    for variant in ['rich', 'chebfix', 'chebest']:
        xf = apply_pc(kept[(variant, 'FS')], r)
        for method in ['MONO-P', 'MONO']:
            xm = apply_pc(kept[(variant, method)], r)
            row = []
            for fn, ix in fields:
                a_ = xm[ix] - xf[ix]; b_ = xf[ix]
                if fn == 'p':
                    a_ = a_ - a_.mean(); b_ = b_ - b_.mean()
                row.append(f"{fn}={np.linalg.norm(a_) / np.linalg.norm(b_):.2e}")
            print(f"   {rname:8s} {variant:8s} {method:7s} " + "  ".join(row))
