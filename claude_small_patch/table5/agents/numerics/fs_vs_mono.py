"""
numerics agent: FS (fieldsplit Schur + MG on HDivSchurPC) vs MONO (MG on the monolithic shifted Jacobian)
on a small mesh, serial.  Standalone copy of the LB_time_slice.py set-up (lines 71-212 of that script,
default flags + the Table-5 flags --length 3e5 --height 4000 --dt 100 --shift 1e-4 --rtol 1e-8 --maxit 100
--reordering [--richardson | --monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left]).
Nothing in the repository is imported except utils.py (the forms), with bytecode writing disabled.

Usage:  python fs_vs_mono.py --nx 4 --nz 4 [--counts_only]
"""
import sys
import os
sys.dont_write_bytecode = True
os.environ.setdefault("OMP_NUM_THREADS", "1")
import json
import re
import time
from argparse import ArgumentParser
from contextlib import ExitStack
from itertools import chain

parser = ArgumentParser()
parser.add_argument('--nx', type=int, default=4)
parser.add_argument('--nz', type=int, default=4)
parser.add_argument('--counts_only', action='store_true', help='only part (c): outer iteration counts')
parser.add_argument('--out', type=str, default=None)
parser.add_argument('--variants', type=str, default='rich,chebfix,chebest,exact,gmres')
parser.add_argument('--methods', type=str, default='FS,FS-P,MONO,MONO-P')
parser.add_argument('--sides', type=str, default='left,right')
args, _ = parser.parse_known_args()

from firedrake import *                      # noqa
from firedrake import dmhooks
import firedrake.utils as fd_utils
from firedrake.petsc import PETSc
import numpy as np
import scipy.sparse as sps

sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
import utils                                  # noqa: E402  (forms of the repository, unchanged)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = args.out or os.path.join(HERE, f"results_nx{args.nx}_nz{args.nz}{'_counts' if args.counts_only else ''}.json")
if args.out is not None and not os.path.isabs(OUT):
    OUT = os.path.join(HERE, OUT)
VIEWDIR = os.path.join(HERE, f"views_nx{args.nx}_nz{args.nz}")
os.makedirs(VIEWDIR, exist_ok=True)

# ----------------------------------------------------------------------------------------------
# Problem set-up: copied from LB_time_slice.py (defaults: use_rotation=False, U_mean=0, degree 2,
# refinement 2) with the Table-5 values.
# ----------------------------------------------------------------------------------------------
nx, nz = args.nx, args.nz
length, height = 3.0e5, 4000.0
refinement, deg = 2, 2
DT, SHIFT = 100.0, 1.0e-4
dt = Constant(DT)
shift = Constant(SHIFT)
U_MEAN = 0.0
use_rotation = False          # LB_time_slice.py: --no_rotation is store_false -> rotation OFF by default
appctx = {"dt": dt, "shift": shift}


def vector_3D(u, uy):
    return u + uy * utils.j()


class HDivSchurPC(AuxiliaryOperatorPC):           # verbatim copy of LB_time_slice.py:103-125
    _prefix = "helmholtzschurpc_"

    def form(self, pc, v, u):
        prefix = (pc.getOptionsPrefix() or "") + self._prefix
        rotation = PETSc.Options().getBool(f"{prefix}use_rotation", False)
        appctx_PC = self.get_appctx(pc)
        dtc = appctx_PC["dt"]
        delta = appctx_PC["shift"]
        n = appctx_PC["n"]
        One = as_vector([1., 1., 1.])
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        unph = Constant(0.5) * (velo + One)
        bnph = Constant(0.5) * (b + Constant(1.))
        w = vector_3D(wxz, wy)
        pnp1 = - Constant(1.) / delta * div(velo)
        Jaux = lhs(utils.LB_velocity(velo, One, unph, w, bnph, pnp1, n, dtc, use_rotation=rotation, twoD=False, U_mean=Constant(U_MEAN)))
        Jaux += lhs(utils.LB_buoyancy(b, Constant(1.), q, unph, bnph, n, dtc, twoD=False, U_mean=Constant(U_MEAN)))
        _, bcs_ = super().form(pc, u, v)
        return (Jaux, bcs_)


class MonolithicShiftPC(AuxiliaryOperatorPC):     # verbatim copy of LB_time_slice.py:127-144
    _prefix = "monoshiftpc_"

    def form(self, pc, v, u):
        test, trial = Jp.arguments()
        a = replace(Jp, {test: v, trial: u})
        _, bcs_ = super().form(pc, u, v)
        return (a, bcs_)


distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (refinement + 1), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]
mesh.name = "finest"

x, y, z = SpatialCoordinate(mesh)
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg - 1)
Pressure = FunctionSpace(mesh, 'DG', deg - 1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure

Un = Function(W)
Unp1 = Function(W)
unxz, uny, bn, pn = split(Un)
unp1xz, unp1y, bnp1, pnp1 = split(Unp1)
w_xz, wy, q, phi = TestFunctions(W)
xc = Constant(length / 2)
a_ic = Constant(5000)
u0_slice, u0yic, b0ic, p0ic = Un.subfunctions
u1_slice, u1yic, b1ic, p1ic = Unp1.subfunctions
b0ic.project(3e-4 * sin(pi * z / height) / (1 + ((x - xc)**2) / a_ic**2))
b1ic.project(3e-4 * sin(pi * z / height) / (1 + ((x - xc)**2) / a_ic**2))
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]
un = vector_3D(unxz, uny)
unp1 = vector_3D(unp1xz, unp1y)
unph = 0.5 * (un + unp1)
bnph = 0.5 * (bn + bnp1)
w = vector_3D(w_xz, wy)
n = FacetNormal(mesh)
appctx.update({"n": n})
eqn = utils.LB_velocity(unp1, un, unph, w, bnph, pnp1, n, dt, use_rotation=use_rotation, U_mean=Constant(U_MEAN))
eqn += utils.LB_buoyancy(bnp1, bn, q, unph, bnph, n, dt, U_mean=Constant(U_MEAN))
eqn += utils.LB_pressure(unp1, phi)
shift_eqn = eqn + shift * pnp1 * phi * dx
Jp = derivative(shift_eqn, Unp1)
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])
U0 = Function(W)
U0.assign(Un)

# Field index sets (serial: field-major numbering)
fis = [iset.indices.copy() for iset in W.dof_dset.field_ises]
IV = np.concatenate(fis[:3])     # (u_xz, u_y, b)  = fieldsplit_1
IP = fis[3]                      # pressure        = fieldsplit_0
NTOT = W.dim()
BCDOFS = np.unique(np.concatenate([bc.nodes for bc in bcs]))   # nodes of W.sub(0) = global ids (field 0 first)
assert np.all(fis[0] == np.arange(len(fis[0])))
print(f"mesh: base nx={nx}, fine {4*nx} x {nz} cells, dofs: u_xz {len(fis[0])}, u_y {len(fis[1])}, b {len(fis[2])}, p {len(fis[3])}, total {NTOT}", flush=True)

# ----------------------------------------------------------------------------------------------
# Solver parameters (Table-5 configuration, LB_time_slice.py:214-398), smoother on MG levels varied.
# ----------------------------------------------------------------------------------------------
PATCH = {   # LB_time_slice.py:220-228 + --reordering (276-280 / 354-358)
    "pc_type": "python",
    "pc_python_type": "firedrake.ASMStarPC",
    "pc_star_construct_dim": 0,
    "pc_star_sub_sub_pc_type": "lu",
    "pc_star_sub_sub_pc_factor_mat_ordering_type": "rcm",
    "pc_star_sub_sub_pc_factor_reuse_ordering": None,
}
SMOOTHERS = {
    'rich':    {'ksp_type': 'richardson', 'ksp_richardson_scale': 0.5, 'ksp_max_it': 6},
    'chebfix': {'ksp_type': 'chebyshev', 'ksp_max_it': 6, 'ksp_chebyshev_eigenvalues': '0.2,2.2'},
    'chebest': {'ksp_type': 'chebyshev', 'ksp_max_it': 6},          # the Table-5 setting
    'gmres':   {'ksp_type': 'gmres', 'ksp_max_it': 6},
}
VARIANT_LABEL = {'rich': '(i) Richardson 0.5 x6', 'chebfix': '(ii) Chebyshev x6 [0.2,2.2] fixed',
                 'chebest': '(iii) Chebyshev x6 PETSc estimate (Table 5)', 'exact': '(iv) exact (LU)',
                 'gmres': '(v) GMRES(6) smoother'}


def mg_block(variant):
    return {'pc_type': 'mg', 'pc_mg_type': 'full', 'pc_mg_cycle_type': 'v',
            'mg_levels': {**PATCH, **SMOOTHERS[variant]},
            'mg_coarse': {'ksp_type': 'preonly', 'pc_type': 'lu'}}


def outer_block(side, outer):
    p = {'snes_type': 'ksponly', 'ksp_type': outer, 'ksp_atol': 0, 'ksp_rtol': 1.0e-8, 'ksp_max_it': 100,
         'ksp_converged_maxits': None, 'ksp_pc_side': side}
    return p


def fs_params(variant, side='left', outer='gmres', use_amat=None):
    p = outer_block(side, outer)
    p.update({'pc_type': 'fieldsplit', 'pc_fieldsplit_type': 'schur', 'pc_fieldsplit_schur_fact_type': 'full',
              'pc_fieldsplit_0_fields': '3', 'pc_fieldsplit_1_fields': '0,1,2',
              'fieldsplit_0': {'ksp_type': 'preonly', 'pc_type': 'bjacobi', 'sub_pc_type': 'ilu'},
              'fieldsplit_1': {'use_rotation': use_rotation, 'ksp_type': 'preonly', 'pc_type': 'python',
                               'pc_python_type': __name__ + '.HDivSchurPC'}})
    if variant == 'exact':
        inner = {'pc_type': 'lu'}
    else:
        inner = mg_block(variant)
        if use_amat is not None:
            inner['pc_use_amat'] = use_amat
    p['fieldsplit_1']['helmholtzschurpc'] = inner
    return p


def mono_params(variant, side='left', outer='gmres', use_amat=None):
    p = outer_block(side, outer)
    if variant == 'exact':
        inner = {'mat_type': 'aij', 'pc_type': 'lu'}
    else:
        inner = {'mat_type': 'aij', **mg_block(variant)}
        if use_amat is not None:
            inner['pc_use_amat'] = use_amat
    p.update({'pc_type': 'python', 'pc_python_type': __name__ + '.MonolithicShiftPC', 'monoshiftpc': inner})
    return p


# method name -> (builder, use_amat).  None = PETSc default (PCMG: use_amat TRUE) = Table-5 configuration.
METHODS = {
    'FS': (fs_params, None),
    'FS-P': (fs_params, False),        # FS with -helmholtzschurpc_pc_use_amat false (MG purely on assembled S_aux)
    'MONO': (mono_params, None),       # Table-5 configuration
    'MONO-P': (mono_params, False),    # MONO with -monoshiftpc_pc_use_amat false (MG purely on A_delta)
}


def make_solver(method, variant, side='left', outer='gmres'):
    builder, ua = METHODS[method]
    params = builder(variant, side=side, outer=outer, use_amat=ua)
    prob = NonlinearVariationalProblem(eqn, Unp1, bcs=bcs, Jp=Jp)
    return NonlinearVariationalSolver(prob, nullspace=nullspace, solver_parameters=params, appctx=appctx)


def solver_context(solver):
    """Replicates the context NonlinearVariationalSolver.solve enters (variational_solver.py:385-397)."""
    stack = ExitStack()
    problem = solver._problem
    forms = (problem.F, problem.J, problem.Jp)
    coefficients = fd_utils.unique(chain.from_iterable(f.coefficients() for f in forms if f is not None))
    solution_dm = solver.snes.getDM()
    problem_dms = []
    for c in chain(coefficients, problem.dirichlet_bcs()):
        dm_ = c.function_space().dm
        if dm_ == solution_dm:
            continue
        if dm_ not in problem_dms:
            problem_dms.append(dm_)
    problem_dms.append(solution_dm)
    for ctx in chain([solver.inserted_options()],
                     [dmhooks.add_hooks(dm_, solver, appctx=solver._ctx) for dm_ in problem_dms],
                     solver._transfer_operators):
        stack.enter_context(ctx)
    return stack


def run_solve(solver):
    Unp1.assign(U0)
    ksp = solver.snes.ksp
    ksp.setConvergenceHistory(length=1000, reset=True)
    t0 = time.perf_counter()
    solver.solve()
    t1 = time.perf_counter()
    its = ksp.getIterationNumber()
    reason = ksp.getConvergedReason()
    hist = [float(v) for v in ksp.getConvergenceHistory()]
    b = ksp.getRhs()
    xs = ksp.getSolution()
    A = ksp.getOperators()[0]
    r = b.duplicate()
    A.mult(xs, r)
    r.aypx(-1.0, b)
    true_rel = r.norm() / b.norm()
    Unp1.assign(U0)
    return dict(its=its, reason=int(reason), hist=hist, true_rel_res=float(true_rel), time=t1 - t0)


def to_dense(M):
    if M.getType() in ('seqaij', 'aij'):
        ai, aj, av = M.getValuesCSR()
        return sps.csr_matrix((av, aj, ai), shape=M.getSize()).toarray()
    mrow, ncol = M.getSize()
    xv, yv = M.createVecs()
    D = np.zeros((mrow, ncol))
    for j in range(ncol):
        xv.set(0.0)
        xv.setValue(j, 1.0)
        xv.assemble()
        M.mult(xv, yv)
        D[:, j] = yv.getArray()
    return D


def ksp_view_text(ksp, name):
    path = os.path.join(VIEWDIR, name + '.txt')
    vw = PETSc.Viewer().createASCII(path, comm=PETSc.COMM_SELF)
    ksp.view(vw)
    vw.destroy()
    with open(path) as f:
        return f.read()


def cheb_info(ksp, name):
    txt = ksp_view_text(ksp, name)
    out = {'type': ksp.getType()}
    mt = re.search(r"eigenvalue targets used: min (\S+), max (\S+)", txt)
    me = re.search(r"eigenvalues estimated via (\S+): min (\S+), max (\S+)", txt)
    if mt:
        out['target'] = [float(mt.group(1).rstrip(',')), float(mt.group(2))]
    if me:
        out['estimated'] = [float(me.group(2).rstrip(',')), float(me.group(3))]
        out['est_ksp'] = me.group(1).rstrip(':')
    return out


def inner_mg(solver, method):
    """Return (inner PETSc PC that is MG or LU, aux python context)."""
    pc = solver.snes.ksp.getPC()
    if method.startswith('FS'):
        ksps = pc.getFieldSplitSubKSP()
        ctx = ksps[1].getPC().getPythonContext()
    else:
        ctx = pc.getPythonContext()
    return ctx.pc, ctx


results = {'nx': nx, 'nz': nz, 'fine_cells': [4 * nx, nz], 'dofs': {'u_xz': len(fis[0]), 'u_y': len(fis[1]),
           'b': len(fis[2]), 'p': len(fis[3]), 'total': NTOT}, 'dt': DT, 'shift': SHIFT}

# ----------------------------------------------------------------------------------------------
# (c) outer iteration counts (and solver objects kept for (a)/(b))
# ----------------------------------------------------------------------------------------------
COUNTS = []
kept = {}
plan = []
VSEL = args.variants.split(',')
MSEL = args.methods.split(',')
SSEL = args.sides.split(',')
for variant in ['rich', 'chebfix', 'chebest', 'exact']:
    for method in ['FS', 'FS-P', 'MONO', 'MONO-P']:
        if variant == 'exact' and method in ('FS-P', 'MONO-P'):
            continue
        for side in ['left', 'right']:
            if variant in VSEL and method in MSEL and side in SSEL:
                plan.append((variant, method, side, 'gmres'))
if 'gmres' in VSEL:
    for method in ['FS', 'FS-P', 'MONO', 'MONO-P']:
        if method in MSEL:
            plan.append(('gmres', method, 'right', 'fgmres'))

T0 = time.perf_counter()
for (variant, method, side, outer) in plan:
    s = make_solver(method, variant, side=side, outer=outer)
    res = run_solve(s)
    res.update(variant=variant, method=method, side=side, outer=outer)
    COUNTS.append(res)
    print(f"[c] {VARIANT_LABEL[variant]:45s} {method:7s} {outer}/{side:5s}: its={res['its']:3d} reason={res['reason']:2d} "
          f"true_rel_res={res['true_rel_res']:.2e}  ({res['time']:.1f}s)", flush=True)
    if side == 'left' or outer == 'fgmres':
        kept[(variant, method)] = s
results['counts'] = [{k: v for k, v in r.items()} for r in COUNTS]
print(f"(c) done in {time.perf_counter() - T0:.1f}s", flush=True)

if args.counts_only:
    with open(OUT, 'w') as f:
        json.dump(results, f, indent=1)
    sys.exit(0)

# ----------------------------------------------------------------------------------------------
# Operators actually used on the MG levels (finest-level Amat / Pmat), Chebyshev intervals
# ----------------------------------------------------------------------------------------------
ref_solver = kept[('chebest', 'FS')]
Jmat, Jpmat = ref_solver.snes.ksp.getOperators()
J_d = to_dense(Jmat)
Jp_d = to_dense(Jpmat)
Dd = Jp_d - J_d
results['norm_Jp_minus_J_pp_over_norm_Jp_pp'] = float(np.linalg.norm(Dd[np.ix_(IP, IP)]) / np.linalg.norm(Jp_d[np.ix_(IP, IP)]))
results['norm_J_pp'] = float(np.linalg.norm(J_d[np.ix_(IP, IP)]))
Dd[np.ix_(IP, IP)] = 0.0
results['norm_Jp_minus_J_outside_pp'] = float(np.linalg.norm(Dd))

ops = {}
intervals = {}
for variant in ['rich', 'chebfix', 'chebest', 'gmres']:
    for method in ['FS', 'FS-P', 'MONO', 'MONO-P']:
        s = kept[(variant, method)]
        mg, ctx = inner_mg(s, method)
        nlev = mg.getMGLevels()
        outerA, outerP = s.snes.ksp.getOperators()
        lev = []
        for l in range(1, nlev):
            kl = mg.getMGSmoother(l)
            A_l, P_l = kl.getOperators()
            info = cheb_info(kl, f"{method}_{variant}_level{l}")
            info.update(level=l, Amat_type=A_l.getType(), Pmat_type=P_l.getType(),
                        Amat_is_outer_J=bool(A_l.handle == outerA.handle),
                        Amat_is_aux=bool(A_l.handle == ctx.P.petscmat.handle),
                        Pmat_is_aux=bool(P_l.handle == ctx.P.petscmat.handle),
                        Amat_has_nullspace=bool(A_l.getNullSpace().handle != 0))
            lev.append(info)
        intervals[f"{method}/{variant}"] = lev
results['levels'] = intervals

# ----------------------------------------------------------------------------------------------
# (a) aux operator == Schur complement
# ----------------------------------------------------------------------------------------------
A00 = Jp_d[np.ix_(IP, IP)]
A01 = Jp_d[np.ix_(IP, IV)]
A10 = Jp_d[np.ix_(IV, IP)]
A11 = Jp_d[np.ix_(IV, IV)]
S_dense = A11 - A10 @ np.linalg.solve(A00, A01)
mg_fs, ctx_fs = inner_mg(kept[('chebest', 'FS')], 'FS')
Saux_d = to_dense(ctx_fs.P.petscmat)
a_res = {'fine_aux_vs_dense_schur_relF': float(np.linalg.norm(Saux_d - S_dense) / np.linalg.norm(S_dense)),
         'fine_norm_S_F': float(np.linalg.norm(S_dense)),
         'fine_norm_graddiv_part_F': float(np.linalg.norm(A10 @ np.linalg.solve(A00, A01))),
         'fine_norm_A11_F': float(np.linalg.norm(A11))}
# MatSchurComplement used as Amat by the FS MG (finest level)
Sschur = mg_fs.getMGSmoother(mg_fs.getMGLevels() - 1).getOperators()[0]
a_res['finest_Amat_type_FS'] = Sschur.getType()
Sschur_d = to_dense(Sschur)
a_res['fine_matschur_vs_dense_schur_relF'] = float(np.linalg.norm(Sschur_d - S_dense) / np.linalg.norm(S_dense))
# level-by-level: Schur complement of MONO-P level operator vs FS level operator
mg_mp, ctx_mp = inner_mg(kept[('chebest', 'MONO-P')], 'MONO-P')
lvl = []
for l in range(mg_fs.getMGLevels()):
    kf = mg_fs.getMGSmoother(l) if l > 0 else mg_fs.getMGCoarseSolve()
    km = mg_mp.getMGSmoother(l) if l > 0 else mg_mp.getMGCoarseSolve()
    Sf = to_dense(kf.getOperators()[1])
    Am = to_dense(km.getOperators()[1])
    nv = Sf.shape[0]
    npl = Am.shape[0] - nv
    a00 = Am[nv:, nv:]
    a01 = Am[nv:, :nv]
    a10 = Am[:nv, nv:]
    a11 = Am[:nv, :nv]
    Sl = a11 - a10 @ np.linalg.solve(a00, a01)
    lvl.append({'level': l, 'n_v': nv, 'n_p': npl, 'relF': float(np.linalg.norm(Sf - Sl) / np.linalg.norm(Sl))})
a_res['levels'] = lvl
results['a'] = a_res
print("(a)", json.dumps(a_res, indent=1), flush=True)

# ----------------------------------------------------------------------------------------------
# (b) preconditioner applications
# ----------------------------------------------------------------------------------------------
rng = np.random.default_rng(20261019)
r_rand = rng.standard_normal(NTOT)
r_rand[BCDOFS] = 0.0
# actual first-step residual F(U^0) (the KSP right-hand side of the first time step)
Fc = assemble(eqn)
with Fc.dat.vec_ro as fv:
    r_res = fv.getArray().copy()
r_res[BCDOFS] = 0.0
# cross-check against the RHS the KSP actually saw
rhs_ksp = None
s_chk = kept[('chebest', 'FS')]
Unp1.assign(U0)
s_chk.snes.ksp.setTolerances(max_it=1)
s_chk.solve()
rhs_ksp = s_chk.snes.ksp.getRhs().getArray().copy()
s_chk.snes.ksp.setTolerances(max_it=100)
Unp1.assign(U0)
results['r_res_vs_ksp_rhs'] = float(np.linalg.norm(rhs_ksp - r_res) / np.linalg.norm(r_res))
results['r_res_norms'] = {'v': float(np.linalg.norm(r_res[IV])), 'p': float(np.linalg.norm(r_res[IP]))}

xs_rand = np.linalg.solve(Jp_d, r_rand)          # x* = A_delta^{-1} r  (LAPACK LU)
xs_res = np.linalg.solve(Jp_d, r_res)
xA_res = np.linalg.lstsq(J_d, r_res, rcond=None)[0]   # unshifted solution (pressure up to a constant)
results['unshifted_solution_residual'] = float(np.linalg.norm(J_d @ xA_res - r_res) / np.linalg.norm(r_res))


def p0(v):           # remove the constant pressure mode (the null space of the unshifted operator)
    return v - v.mean()


def apply_pc(solver, r):
    pc = solver.snes.ksp.getPC()
    xv, yv = solver.snes.ksp.getOperators()[0].createVecs()
    xv.setArray(r.copy())
    yv.set(0.0)
    with solver_context(solver):
        pc.apply(xv, yv)
    return yv.getArray().copy()


def stats(xm, r, xref, xA=None):
    xv_, xp_ = xm[IV], xm[IP]
    d = {}
    d['e_v'] = float(np.linalg.norm(xv_ - xref[IV]) / np.linalg.norm(xref[IV]))
    d['e_p'] = float(np.linalg.norm(p0(xp_ - xref[IP])) / np.linalg.norm(p0(xref[IP])))
    backsub = np.linalg.solve(A00, r[IP] - A01 @ xv_)          # (delta M)^{-1} (r_p - D x_v)
    d['z_rel'] = float(np.linalg.norm(p0(xp_ - backsub)) / np.linalg.norm(p0(xp_)))
    d['z_raw_rel'] = float(np.linalg.norm(xp_ - backsub) / np.linalg.norm(xp_))
    d['p_mean'] = float(xp_.mean())
    rs = r - Jp_d @ xm                                         # shifted-system residual after one application
    ru = r - J_d @ xm                                          # unshifted residual after one application
    nr = np.linalg.norm(r)
    d['shift_res_v'] = float(np.linalg.norm(rs[IV]) / nr)
    d['shift_res_p'] = float(np.linalg.norm(rs[IP]) / nr)
    d['unsh_res_v'] = float(np.linalg.norm(ru[IV]) / nr)
    d['unsh_res_p'] = float(np.linalg.norm(ru[IP]) / nr)
    if xA is not None:
        d['eA_v'] = float(np.linalg.norm(xv_ - xA[IV]) / np.linalg.norm(xA[IV]))
        d['eA_p'] = float(np.linalg.norm(p0(xp_ - xA[IP])) / np.linalg.norm(p0(xA[IP])))
    return d


b_res = {}
for rname, r, xref, xA in [('random', r_rand, xs_rand, None), ('residual', r_res, xs_res, xA_res)]:
    b_res[rname] = {}
    # reference rows
    b_res[rname]['xstar_A_vs_xstar_Adelta'] = None if xA is None else {
        'v': float(np.linalg.norm(xA[IV] - xref[IV]) / np.linalg.norm(xref[IV])),
        'p': float(np.linalg.norm(p0(xA[IP] - xref[IP])) / np.linalg.norm(p0(xref[IP])))}
    for variant in ['rich', 'chebfix', 'chebest', 'exact', 'gmres']:
        outs = {}
        for method in ['FS', 'FS-P', 'MONO', 'MONO-P']:
            if (variant, method) in kept:
                outs[method] = apply_pc(kept[(variant, method)], r)
        if variant == 'exact':
            outs['FS-P'] = outs['FS']
            outs['MONO-P'] = outs['MONO']
        rowv = {}
        xf = outs['FS']
        for method, xm in outs.items():
            st = stats(xm, r, xref, xA)
            st['dv_vs_FS'] = float(np.linalg.norm(xm[IV] - xf[IV]) / np.linalg.norm(xf[IV]))
            st['dp_vs_FS'] = float(np.linalg.norm(p0(xm[IP] - xf[IP])) / np.linalg.norm(p0(xf[IP])))
            rowv[method] = st
        # where the pressure difference lives: x_p^FS - x_p^M  vs  -(z of M)
        for method in ['MONO', 'MONO-P']:
            xm = outs[method]
            backsub = np.linalg.solve(A00, r[IP] - A01 @ xm[IV])
            zM = p0(xm[IP] - backsub)
            dpf = p0(xf[IP] - xm[IP])
            rowv[method]['dp_minus_minus_z_rel'] = float(np.linalg.norm(dpf + zM) / max(np.linalg.norm(dpf), 1e-300))
        b_res[rname][variant] = rowv
        print(f"(b) {rname} {VARIANT_LABEL[variant]}", flush=True)
        for method, st in rowv.items():
            print(f"     {method:7s} dv={st['dv_vs_FS']:.2e} dp={st['dp_vs_FS']:.2e} e_v={st['e_v']:.2e} e_p={st['e_p']:.2e} "
                  f"z={st['z_rel']:.2e} shres(v,p)=({st['shift_res_v']:.2e},{st['shift_res_p']:.2e}) "
                  f"unres(v,p)=({st['unsh_res_v']:.2e},{st['unsh_res_p']:.2e})"
                  + (f" eA(v,p)=({st['eA_v']:.2e},{st['eA_p']:.2e})" if 'eA_v' in st else ''), flush=True)
results['b'] = b_res

# ----------------------------------------------------------------------------------------------
# Mechanism checks on the finest level of MONO-P (pure A_delta MG)
# ----------------------------------------------------------------------------------------------
mech = {}
mg_mpf, ctx_mpf = inner_mg(kept[('chebfix', 'MONO-P')], 'MONO-P')
kfin = mg_mpf.getMGSmoother(mg_mpf.getMGLevels() - 1)
asm_pc = kfin.getPC()
Afin = kfin.getOperators()[1]


def asm_apply(v):
    xv_, yv_ = Afin.createVecs()
    xv_.setArray(v.copy())
    with solver_context(kept[('chebfix', 'MONO-P')]):
        asm_pc.apply(xv_, yv_)
    return yv_.getArray().copy()


# (M1) pressure rows of A_delta * ASM^{-1}:  (A_delta M^{-1} r)_p == 2 r_p
y = asm_apply(r_rand)
AMr = Jp_d @ y
mech['M1_AdeltaMinv_pp_minus_2I_rel'] = float(np.linalg.norm(AMr[IP] - 2.0 * r_rand[IP]) / np.linalg.norm(2.0 * r_rand[IP]))
# velocity of ASM output depends on r only through the reduced residual  r_v - A10 A00^{-1} r_p
r2 = r_rand.copy()
r2[IP] = 0.0
r2[IV] = r_rand[IV] - A10 @ np.linalg.solve(A00, r_rand[IP])
y2 = asm_apply(r2)
mech['M1b_asm_velocity_depends_on_reduced_residual_only_rel'] = float(np.linalg.norm(y2[IV] - y[IV]) / np.linalg.norm(y[IV]))


# (M2) Chebyshev [0.2,2.2] x6 (PETSc first-kind recurrence, zero guess): pressure residual = q6(2) r_p
def cheb_numpy(Aop, Minv, b, emin, emax, its):
    scale = 2.0 / (emax + emin)
    alpha = 1.0 - scale * emin
    mu = 1.0 / alpha
    omegaprod = 2.0 / alpha
    c_km1, c_k = 1.0, mu
    p_km1 = np.zeros_like(b)
    p_k = scale * Minv(b) + p_km1
    for i in range(1, its):
        rr = b - Aop(p_k)
        zz = Minv(rr)
        c_kp1 = 2.0 * mu * c_k - c_km1
        omega = omegaprod * c_k / c_kp1
        p_kp1 = (1.0 - omega) * p_km1 + omega * p_k + omega * scale * zz
        p_km1, p_k = p_k, p_kp1
        c_km1, c_k = c_k, c_kp1
    return p_k


def T6(t):
    return 32 * t**6 - 48 * t**4 + 18 * t**2 - 1


q6_2 = T6((2.2 + 0.2 - 2 * 2.0) / (2.2 - 0.2)) / T6((2.2 + 0.2) / (2.2 - 0.2))
xc6 = cheb_numpy(lambda v: Jp_d @ v, asm_apply, r_rand, 0.2, 2.2, 6)
rp_after = (r_rand - Jp_d @ xc6)[IP]
mech['M2_q6_at_2_formula'] = float(q6_2)
mech['M2_numpy_cheb_pressure_residual_ratio'] = float(np.dot(rp_after, r_rand[IP]) / np.dot(r_rand[IP], r_rand[IP]))
mech['M2_numpy_cheb_pressure_residual_minus_q6r_rel'] = float(np.linalg.norm(rp_after - q6_2 * r_rand[IP]) / np.linalg.norm(q6_2 * r_rand[IP]))
# the same with the PETSc smoother object itself (has the constant-pressure null space attached -> constant removed)
bv, xv0 = Afin.createVecs()[1], Afin.createVecs()[0]
bv.setArray(r_rand.copy())
xv0.set(0.0)
with solver_context(kept[('chebfix', 'MONO-P')]):
    kfin.solve(bv, xv0)
xp6 = xv0.getArray().copy()
rp_after2 = (r_rand - Jp_d @ xp6)[IP]
mech['M2_petsc_smoother_vs_numpy_cheb_rel'] = float(np.linalg.norm(p0(xp6[IP] - xc6[IP])) / np.linalg.norm(p0(xc6[IP])))
mech['M2_petsc_smoother_vs_numpy_cheb_vel_rel'] = float(np.linalg.norm(xp6[IV] - xc6[IV]) / np.linalg.norm(xc6[IV]))
# Richardson 0.5 one step: pressure residual exactly zero
x1 = 0.5 * asm_apply(r_rand)
mech['M2_rich05_one_step_pressure_residual_rel'] = float(np.linalg.norm((r_rand - Jp_d @ x1)[IP]) / np.linalg.norm(r_rand[IP]))
results['mechanism'] = mech
print("(mech)", json.dumps(mech, indent=1), flush=True)

with open(OUT, 'w') as f:
    json.dump(results, f, indent=1)
print("wrote", OUT, flush=True)
