"""Build lb_audit_copy.py = LB_time_slice.py (unchanged repo file, read only) + the proposed Q7 flags + an audit dump.
Every insertion is anchored on an exact line of the original and asserted to occur exactly once."""
import pathlib, sys
SRC = pathlib.Path('/home/mh2520/molin/phd/yr2/Boussinesq_implicit/LB_time_slice.py')
DST = pathlib.Path(__file__).with_name('lb_audit_copy.py')
s = SRC.read_text()

def insert_after(s, anchor, text):
    assert s.count(anchor) == 1, anchor
    return s.replace(anchor, anchor + text)

# 0. find utils.py of the repo
s = 'import sys; sys.path.insert(0, "/home/mh2520/molin/phd/yr2/Boussinesq_implicit")\n' + s

# 1. flags (after the last existing add_argument)
FLAGS = r'''
parser.add_argument('--inner_lu', action='store_true', help='Diagnostic: exact LU inner solves in both solvers (FS: LU of the assembled Schur auxiliary operator S_delta instead of the MG cycle; MONO: LU of A_delta instead of the MG cycle).')
parser.add_argument('--fs_smoother', type=str, default=None, choices=['chebyshev', 'richardson', 'gmres'], help='Diagnostic: KSP around the patch smoother on the MG levels of the fieldsplit solver, with fieldsplit_1 kept preonly (ONE MG cycle per outer iteration). Default: chebyshev with --richardson, gmres (and fieldsplit_1 fgmres) without.')
parser.add_argument('--richardson_scale', type=float, default=None, help='Diagnostic: damping of a Richardson MG level smoother, both solvers (PETSc default 1.0).')
parser.add_argument('--cheb_interval', type=float, nargs=2, default=None, metavar=('EMIN', 'EMAX'), help='Diagnostic: fixed Chebyshev interval on every MG level of both solvers (-ksp_chebyshev_eigenvalues; switches the GMRES eigenvalue estimate off).')
parser.add_argument('--mg_use_pmat', action='store_true', help='Diagnostic: finest-level smoother and residual of the MG cycle on the auxiliary (shifted) operator instead of the outer Amat (PCMG default useAmat=TRUE: unshifted J for MONO, the shifted Schur complement for FS).')
'''
s = insert_after(s, "parser.add_argument('--outer_ksp', type=str, default=None, choices=['gmres', 'fgmres'], help='Diagnostic: outer Krylov method (the script\\'s defaults: gmres for the fieldsplit solver with --richardson, fgmres otherwise and for the monolithic solver; the 2026-06-01 monolithic runs used gmres with left preconditioning).')\n", FLAGS)

# 2. FS inner-solver overrides, before params_schur is built
FS_OVR = r'''
if args.inner_lu:   # (i) exact LU of the assembled S_delta instead of the MG cycle
    helmholtz_schur_pc_params = {'pc_type': 'lu', 'pc_factor_mat_solver_type': 'mumps'}
elif not args.direct:
    if args.fs_smoother is not None:          # (ii)/(iv) level KSP of the FS cycle
        helmholtz_schur_pc_params['mg_levels_ksp_type'] = args.fs_smoother
    if args.richardson_scale is not None:     # (ii)
        helmholtz_schur_pc_params['mg_levels_ksp_richardson_scale'] = args.richardson_scale
    if args.cheb_interval is not None:        # (iii)
        helmholtz_schur_pc_params['mg_levels_ksp_chebyshev_eigenvalues'] = '%r,%r' % tuple(args.cheb_interval)
    if args.mg_use_pmat:                      # (v)
        helmholtz_schur_pc_params['pc_use_amat'] = False
'''
anchor = "\nparams_schur = {\n"
assert s.count(anchor) == 1
s = s.replace(anchor, FS_OVR + anchor)

# 3. keep fieldsplit_1 preonly when a level smoother is chosen explicitly (before --inner_rtol, which wins)
s = insert_after(s, "params_schur['fieldsplit_1_ksp_type'] = 'preonly' if args.richardson else 'fgmres'\n",
                 "if args.fs_smoother is not None:\n    params_schur['fieldsplit_1_ksp_type'] = 'preonly'\n")

# 4. MONO level smoother overrides
MONO_LEV = r'''
if args.richardson_scale is not None:         # (ii)
    mg_levels_mono['ksp_richardson_scale'] = args.richardson_scale
if args.cheb_interval is not None:            # (iii)
    mg_levels_mono['ksp_chebyshev_eigenvalues'] = '%r,%r' % tuple(args.cheb_interval)
'''
s = insert_after(s, "        'pc_star_sub_sub_pc_factor_reuse_ordering': None,\n    })\n\nparams_monolithic = {", "")  # sanity only
s = s.replace("        'pc_star_sub_sub_pc_factor_reuse_ordering': None,\n    })\n\nparams_monolithic = {",
              "        'pc_star_sub_sub_pc_factor_reuse_ordering': None,\n    })\n" + MONO_LEV + "\nparams_monolithic = {")

# 5. MONO inner overrides (after the --inner_rtol block)
MONO_OVR = r'''if args.inner_lu:   # (i) exact LU of A_delta instead of the MG cycle
    params_monolithic['monoshiftpc'] = {'mat_type': 'aij', 'pc_type': 'lu', 'pc_factor_mat_solver_type': 'mumps'}
elif args.mg_use_pmat and args.inner_rtol is None:   # (v)
    params_monolithic['monoshiftpc']['pc_use_amat'] = False
'''
anchor = "if args.pc_side is not None:\n    params_monolithic['ksp_pc_side'] = args.pc_side\n"
assert s.count(anchor) == 1
s = s.replace(anchor, MONO_OVR + anchor)

# 6. audit dump appended at the end (only with AUDIT_DUMP=<file.npz>)
s += r'''

# ============================ code_audit: dump operators =============================
import os
if os.environ.get("AUDIT_DUMP"):
    import numpy as np
    def dense(M):   # NB: petsc4py Mat.convert() without 'out' converts IN PLACE; read the values instead
        m, n = M.getSize()
        return M.getValues(range(m), range(n))
    def python_dense(M):
        m, n = M.getSize(); x, y = M.createVecs()
        D = np.zeros((m, n))
        for j in range(n):
            x.set(0.0); x.setValue(j, 1.0); x.assemble(); M.mult(x, y); D[:, j] = y.getArray()
        return D
    def mg_dump(mg, tag, out):
        nlev = mg.getMGLevels()
        for l in range(nlev):
            k = mg.getMGCoarseSolve() if l == 0 else mg.getMGSmoother(l)
            A, P = k.getOperators()
            out[f'{tag}_L{l}_P'] = dense(P)
            out[f'{tag}_L{l}_Atype'] = np.array(A.getType()); out[f'{tag}_L{l}_A_is_P'] = np.array(A.handle == P.handle)
            out[f'{tag}_L{l}_A_is_J'] = np.array(A.handle == nsolver._ctx._jac.petscmat.handle)
            out[f'{tag}_L{l}_A_null'] = np.array(A.getNullSpace().handle != 0); out[f'{tag}_L{l}_P_null'] = np.array(P.getNullSpace().handle != 0)
            if l > 0:
                out[f'{tag}_L{l}_I'] = python_dense(mg.getMGInterpolation(l))
                asm = k.getPC().getPythonContext().asmpc
                subs = asm.getASMSubKSP()
                pm = [dense(sk.getOperators()[0]) for sk in subs]
                out[f'{tag}_L{l}_npatch'] = np.array(len(pm)); out[f'{tag}_L{l}_patchsize'] = np.array([p.shape[0] for p in pm])
                # exactness of the patch LU (preonly + lu, rcm) against a dense solve
                rng = np.random.default_rng(0); errs = []; conds = []
                for sk, p in zip(subs, pm):
                    b = sk.getOperators()[0].createVecLeft(); xk = b.duplicate()
                    bb = rng.standard_normal(p.shape[0]); b.setArray(bb); sk.solve(b, xk)
                    xr = np.linalg.solve(p, bb); errs.append(np.linalg.norm(xk.getArray() - xr) / np.linalg.norm(xr)); conds.append(np.linalg.cond(p))
                out[f'{tag}_L{l}_patch_lu_relerr'] = np.array(errs); out[f'{tag}_L{l}_patch_cond'] = np.array(conds)
                out[f'{tag}_L{l}_patches'] = np.array(pm)
    out = {}
    J = nsolver._ctx._jac.petscmat; Jp = nsolver._ctx._pjac.petscmat
    out['J'] = dense(J); out['Jp'] = dense(Jp)
    out['J_null'] = np.array(J.getNullSpace().handle != 0); out['Jp_null'] = np.array(Jp.getNullSpace().handle != 0)
    out['ises'] = np.array([len(i.indices) for i in W._ises]); out['ises_first'] = np.array([i.indices[0] for i in W._ises])
    out['bc_nodes'] = np.array(sorted(set(bc1.nodes) | set(bc2.nodes)))
    out['use_rotation'] = np.array(use_rotation)
    pc = nsolver.snes.ksp.getPC()
    if args.monolithic:
        ctx = pc.getPythonContext(); out['aux'] = dense(ctx.P.petscmat)
        out['aux_null'] = np.array(ctx.P.petscmat.getNullSpace().handle != 0)
        out['inner_A_is_J'] = np.array(ctx.pc.getOperators()[0].handle == J.handle)
        if ctx.pc.getType() == 'mg':
            mg_dump(ctx.pc, 'mg', out)
    else:
        ksp0, ksp1 = pc.getFieldSplitSubKSP()
        A0, P0 = ksp0.getOperators(); out['fs0_P'] = dense(P0); out['fs0_A_null'] = np.array(A0.getNullSpace().handle != 0)
        rng = np.random.default_rng(1); x = P0.createVecLeft(); y = x.duplicate(); xx = rng.standard_normal(x.getSize()); x.setArray(xx)
        ksp0.solve(x, y); out['fs0_x'] = xx; out['fs0_y'] = y.getArray().copy()
        ctx = ksp1.getPC().getPythonContext(); out['aux'] = dense(ctx.P.petscmat)
        out['aux_null'] = np.array(ctx.P.petscmat.getNullSpace().handle != 0)
        A1, P1 = ksp1.getOperators(); out['fs1_Atype'] = np.array(A1.getType()); out['fs1_P'] = dense(P1)
        out['inner_A_is_S'] = np.array(ctx.pc.getOperators()[0].handle == A1.handle)
        if ctx.pc.getType() == 'mg':
            mg_dump(ctx.pc, 'mg', out)
    np.savez(os.environ["AUDIT_DUMP"], **out)
    print("AUDIT_DUMP written", os.environ["AUDIT_DUMP"])
'''
DST.write_text(s)
print('written', DST)
