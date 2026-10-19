#!/usr/bin/env python
'''
inner_solve.py
==============
Isolated "field 1" solve of Nonlinear_Boussinesq_Irk_SC.py.

The solver script preconditions the outer FGMRES by an auxiliary operator (HDivSchurPC)
followed by PCFIELDSPLIT (Schur, full).  Field 1 = (u, u_y, b) is handled by ONE
application of geometric multigrid.  This script takes that multigrid out of the outer
Schur-complement iteration and measures it directly:

    solve   A11 x = rhs     with FGMRES (rtol 1e-8), preconditioned by exactly the
                            multigrid of the solver script (pc_mg_type full, V cycle,
                            GMRES(6) + ASM patch smoother with LU patch solves (rcm),
                            coarse LU), for the star patch and for the vanka patch.

A11 is the (u, u_y, b) block of the auxiliary operator, linearised about the state of
rest, f = 0.  With the stage derivative k = (k_u, k_uy, k_b) of the implicit midpoint
rule as unknown (state = state_n + dt/2 * k):

    a(k, test) = (w, k_u) + (w_y, k_uy) + (dt / (2 delta)) (div w, div k_u)
                 - (dt/2) (w . zhat, k_b) + (q, k_b) + (dt/2) N^2 (q, zhat . k_u)

Eliminating k_b gives the weighted H(div) form
    (u_x, w_x) + m (u_z, w_z) + c (div u, div w),  m = 1 + (dt N)^2/4,  c = dt/(2 delta).

Mesh hierarchy, function spaces and boundary conditions are built with the same calls
as the solver script (the functions of the repository's utils.py are imported, nothing
in the repository is modified and no byte code is written there).

Usage: every physical / numerical parameter accepts a LIST; the script runs the Cartesian
product of all lists inside ONE python process (meshes and compiled forms are reused).

    python inner_solve.py --tag K --K 0.1 1 10 100 1000 --patch star vanka --rhs random bump

The shift is delta = K * delta_thr (claude_paper_shift.py formula) unless --c is given,
in which case c = dt / (2 delta) is prescribed directly (delta = dt / (2 c)).

Output: one line per solve appended to <out>.csv, residual histories to <out>_hist.jsonl.
'''
import os
import sys
import time
import json
import csv
import argparse
import itertools

sys.dont_write_bytecode = True          # do not write __pycache__ into the repository
REPO = os.environ.get('BOUSSINESQ_REPO', '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
if REPO not in sys.path:
    sys.path.insert(0, REPO)

import numpy as np                                                   # noqa: E402
from firedrake import *                                              # noqa: E402,F403
from firedrake.preconditioners.asm import ASMStarPC, ASMVankaPC     # noqa: E402
from firedrake.petsc import PETSc                                    # noqa: E402
import utils                                                         # noqa: E402  (repository, read only)

Print = PETSc.Sys.Print
N_FREQ = 1.0e-2           # = sqrt(utils.buo_freq())


# ----------------------------------------------------------------------------------
#  Shift formula (copy of claude_paper_shift.py, kept here so that the script is
#  self-contained; checked against the module in verify_against_solver.py)
# ----------------------------------------------------------------------------------
def m_factor(dt, N=N_FREQ):
    return 1.0 + (dt * N) ** 2 / 4.0


def C_M_squared(dt, Lx, Lz, N=N_FREQ):
    return max((Lx / (2.0 * np.pi)) ** 2, m_factor(dt, N) * (Lz / np.pi) ** 2)


def regime(dt, Lx, Lz, N=N_FREQ):
    return 'weak' if m_factor(dt, N) * (Lz / np.pi) ** 2 <= (Lx / (2.0 * np.pi)) ** 2 else 'strong'


def shift_threshold(dt, Lx, Lz, N=N_FREQ):
    return dt / (2.0 * C_M_squared(dt, Lx, Lz, N))


# ----------------------------------------------------------------------------------
#  Patch preconditioners with counters.  They ARE firedrake.ASMStarPC / ASMVankaPC
#  (same options prefixes pc_star_ / pc_vanka_); only the number of applications and
#  the patch sizes are recorded.
# ----------------------------------------------------------------------------------
class _CountingMixin(object):
    def get_patches(self, V):
        ises = super().get_patches(V)
        self._patch_sizes = [int(i.getSize()) for i in ises]
        return ises

    def initialize(self, pc):
        self._napply = 0
        super().initialize(pc)

    def apply(self, pc, x, y):
        self._napply += 1
        super().apply(pc, x, y)


class CountingASMStarPC(_CountingMixin, ASMStarPC):
    pass


class CountingASMVankaPC(_CountingMixin, ASMVankaPC):
    pass


def patch_parameters(patch):
    '''The patch dictionaries of Nonlinear_Boussinesq_Irk_SC.py (only the python type is
    replaced by the counting subclass).'''
    if patch == 'star':
        return {
            'pc_type': 'python',
            'pc_python_type': __name__ + '.CountingASMStarPC',
            'pc_star_construct_dim': 0,
            'pc_star_sub_sub_pc_type': 'lu',
            'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
            'pc_star_sub_sub_pc_factor_reuse_ordering': None,
        }
    elif patch == 'vanka':
        return {
            'pc_type': 'python',
            'pc_python_type': __name__ + '.CountingASMVankaPC',
            'pc_vanka_construct_dim': 1,
            'pc_vanka_sub_sub_pc_type': 'lu',
            'pc_vanka_sub_sub_pc_factor_mat_ordering_type': 'rcm',
            'pc_vanka_sub_sub_pc_factor_reuse_ordering': None,
        }
    raise ValueError(patch)


def mg_parameters(patch, smooth_its=6, smoother_ksp='gmres', mg_type='full'):
    '''The fieldsplit_1 dictionary of Nonlinear_Boussinesq_Irk_SC.py.'''
    return {
        'pc_type': 'mg',
        'pc_mg_type': mg_type,
        'pc_mg_cycle_type': 'v',
        'mg_levels': {
            'ksp_type': smoother_ksp,
            'ksp_max_it': smooth_its,
            **patch_parameters(patch),
        },
        'mg_coarse': {
            'ksp_type': 'preonly',
            'pc_type': 'lu',
        },
    }


def outer_parameters(patch, smooth_its=6, rtol=1.0e-8, maxit=200, restart=200,
                     smoother_ksp='gmres', mg_type='full', outer='fgmres'):
    params = {
        'mat_type': 'aij',
        'snes_type': 'ksponly',
        'ksp_type': outer,
        'ksp_rtol': rtol,
        'ksp_atol': 0.0,
        'ksp_max_it': maxit,
        'ksp_gmres_restart': restart,
        'ksp_converged_maxits': None,        # do not raise when the cap is reached
    }
    if outer == 'fgmres':
        params['ksp_pc_side'] = 'right'
    if outer == 'richardson':
        params['ksp_norm_type'] = 'unpreconditioned'
    params.update(mg_parameters(patch, smooth_its, smoother_ksp, mg_type))
    return params


# ----------------------------------------------------------------------------------
#  Mesh hierarchy and spaces: the lines of Nonlinear_Boussinesq_Irk_SC.py
# ----------------------------------------------------------------------------------
def vector_3D(u, uy):
    return u + uy * utils.j()


def build_hierarchy(nx, nz, Lx, Lz, refinement=2):
    distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
    m = PeriodicIntervalMesh(nx, Lx, distribution_parameters=distribution_parameters)
    mh = MeshHierarchy(m, refinement_levels=refinement)
    hierarchy = ExtrudedMeshHierarchy(mh, Lz, layers=[nz] * (refinement + 1), extrusion_type='uniform')
    new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
    return new_mh


def build_space(mesh, deg=2):
    '''Field 1 of the solver script: W.sub(0) * W.sub(1) * W.sub(2).'''
    V_2D = utils.extrude_RT(mesh, k=deg)
    Vy = FunctionSpace(mesh, 'DG', deg - 1)
    Vb = utils.W_theta(mesh, k=deg)
    return V_2D * Vy * Vb


def boundary_conditions(W1):
    bc1 = DirichletBC(W1.sub(0), as_vector([0., 0., 0.]), "top")
    bc2 = DirichletBC(W1.sub(0), as_vector([0., 0., 0.]), "bottom")
    return [bc1, bc2]


# ----------------------------------------------------------------------------------
#  Bilinear forms
# ----------------------------------------------------------------------------------
def field1_form_rest(W1, dt, delta, N2=None):
    '''(u, u_y, b) block of the auxiliary operator linearised about rest, f = 0.'''
    N2 = utils.buo_freq() if N2 is None else N2
    uxz, uy, b = TrialFunctions(W1)
    wxz, wy, q = TestFunctions(W1)
    u = vector_3D(uxz, uy)
    w = vector_3D(wxz, wy)
    kz = utils.k()
    half = Constant(0.5)
    a = (
        inner(w, u) * dx
        + half * dt / delta * div(w) * div(u) * dx
        - half * dt * inner(w, kz) * b * dx
        + q * b * dx
        + half * dt * N2 * q * inner(kz, u) * dx
    )
    return a


def field1_form_utils(W1, dt_value, delta, state=None, stage=None, use_rotation=False):
    '''
    The same block built the way the solver builds it: HDivSchurPC.getNewForm (with the
    pressure dropped, p -> -(1/delta) div u) followed by Irksome's stage form for
    GaussLegendre(1) and the derivative with respect to the stage.  `state` is the
    linearisation state (u, u_y, b) at the beginning of the step, `stage` the stage
    derivative about which the Jacobian is taken.  Returns (J, stage_bcs, state, stage).
    '''
    from irksome import GaussLegendre, MeshConstant
    from irksome.stage_derivative import getForm
    mesh = W1.mesh()
    MC = MeshConstant(mesh)
    dtc = MC.Constant(dt_value)
    t = MC.Constant(0.0)
    n = FacetNormal(mesh)
    state = Function(W1) if state is None else state
    stage = Function(W1) if stage is None else stage
    test = TestFunction(W1)
    uxz, uy, b = split(state)
    wxz, wy, q = split(test)
    u = vector_3D(uxz, uy)
    w = vector_3D(wxz, wy)
    p = - Constant(1.) / delta * div(u)
    F = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=use_rotation)
    F += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
    bcs = boundary_conditions(W1)
    Fs, bcs_s = getForm(F, GaussLegendre(1), t, dtc, state, stage, bcs=bcs)
    J = derivative(Fs, stage, du=TrialFunction(W1))
    return J, bcs_s, state, stage


def initial_buoyancy(W1, Lx, Lz):
    '''b0 of the solver script.'''
    mesh = W1.mesh()
    x, y, z = SpatialCoordinate(mesh)
    xc = Constant(Lx / 2)
    a = Constant(5000)
    b0 = Function(W1.sub(2).collapse() if False else W1.sub(2))
    b0 = Function(W1.sub(2).collapse())
    b0.project(3e-4 * sin(pi * z / Lz) / (1 + ((x - xc) ** 2) / a ** 2))
    return b0


# ----------------------------------------------------------------------------------
#  Right-hand sides
# ----------------------------------------------------------------------------------
def make_rhs(kind, W1, Lx, Lz, seed=1234):
    mesh = W1.mesh()
    wxz, wy, q = TestFunctions(W1)
    w = vector_3D(wxz, wy)
    x, y, z = SpatialCoordinate(mesh)
    if kind == 'bump':
        # residual of the first Newton step of the first time step of the solver script:
        # state at rest with the buoyancy perturbation b0  ->  rhs = (w . zhat, b0)
        b0 = initial_buoyancy(W1, Lx, Lz)
        L = inner(w, utils.k()) * b0 * dx
        return assemble(L)
    if kind == 'smooth':
        kx = 2 * pi / Lx
        kzz = pi / Lz
        f = as_vector([sin(kx * x) * cos(kzz * z), cos(kx * x) * cos(kzz * z), cos(kx * x) * sin(kzz * z)])
        g = Constant(N_FREQ) * cos(kx * x) * sin(kzz * z)
        L = inner(w, f) * dx + q * g * dx
        return assemble(L)
    if kind == 'random':
        rhs = Cofunction(W1.dual())
        rank = mesh.comm.rank
        rng = np.random.default_rng(seed + 7919 * rank)
        for d in rhs.dat:
            d.data[:] = rng.standard_normal(d.data.shape)
        return rhs
    raise ValueError(kind)


# ----------------------------------------------------------------------------------
#  One solve
# ----------------------------------------------------------------------------------
_solver_counter = [0]


def collect_patch_stats(ksp):
    '''Walk the PCMG of `ksp` and return per-level statistics of the ASM smoother.'''
    pc = ksp.getPC()
    comm = ksp.comm.tompi4py()
    stats = []
    if pc.getType() != 'mg':
        return stats
    nlev = pc.getMGLevels()
    for lev in range(1, nlev):
        sm = pc.getMGSmoother(lev)
        ctx = sm.getPC().getPythonContext()
        A, P = sm.getOperators()
        sizes = getattr(ctx, '_patch_sizes', [])
        nz_fact = 0.0
        nz_mat = 0.0
        nsub = 0
        if ctx.asmpc.getType() == 'asm':
            for sub in ctx.asmpc.getASMSubKSP():
                S, _ = sub.getOperators()
                if S.getSize()[0] == 0:
                    continue
                nsub += 1
                nz_mat += S.getInfo()['nz_used']
                try:
                    nz_fact += sub.getPC().getFactorMatrix().getInfo()['nz_used']
                except Exception:
                    nz_fact += float('nan')
        from mpi4py import MPI
        stats.append(dict(
            level=lev,
            napply=int(ctx._napply),
            npatch=int(comm.allreduce(nsub, op=MPI.SUM)),
            patch_size_max=int(comm.allreduce(max(sizes) if sizes else 0, op=MPI.MAX)),
            patch_size_min=int(comm.allreduce(min(sizes) if sizes else 10**9, op=MPI.MIN)),
            nz_patch_factors=float(comm.allreduce(nz_fact, op=MPI.SUM)),
            nz_patch_mats=float(comm.allreduce(nz_mat, op=MPI.SUM)),
            nz_level_mat=float(P.getInfo(PETSc.Mat.InfoType.GLOBAL_SUM)['nz_used']),
            ndof_level=int(P.getSize()[0]),
        ))
    return stats


def run_case(cache, nx, nz, Lx, Lz, dt, K=None, c=None, patch='star', smooth_its=6, rhs='random',
             state='rest', refinement=2, degree=2, rtol=1.0e-8, maxit=200, restart=200,
             smoother_ksp='gmres', mg_type='full', outer='fgmres', timing_repeats=1,
             view_file=None, tag=''):
    comm = COMM_WORLD
    # ---- shift ----
    m = m_factor(dt)
    CM2 = C_M_squared(dt, Lx, Lz)
    if c is not None:
        delta = dt / (2.0 * c)
        K_eff = CM2 / c
    else:
        delta = K * shift_threshold(dt, Lx, Lz)
        K_eff = K
        c = dt / (2.0 * delta)
    ncol = nx * 2 ** refinement
    dxf = Lx / ncol
    dzf = Lz / nz

    # ---- mesh, space (cached) ----
    key = (nx, nz, Lx, Lz, refinement, degree)
    if key not in cache:
        cache.clear()                     # keep at most one hierarchy alive
        mh = build_hierarchy(nx, nz, Lx, Lz, refinement)
        W1 = build_space(mh[-1], degree)
        cache[key] = dict(mh=mh, W1=W1, rhs={})
    W1 = cache[key]['W1']
    if rhs not in cache[key]['rhs']:
        cache[key]['rhs'][rhs] = make_rhs(rhs, W1, Lx, Lz)
    L = cache[key]['rhs'][rhs]

    dtc = Constant(dt)
    deltac = Constant(delta)
    if state == 'rest':
        a = field1_form_rest(W1, dtc, deltac)
        bcs = boundary_conditions(W1)
    elif state == 'ic':
        st = Function(W1)
        st.subfunctions[2].assign(initial_buoyancy(W1, Lx, Lz))
        a, bcs, _, _ = field1_form_utils(W1, dt, deltac, state=st)
    elif state == 'rest_utils':
        a, bcs, _, _ = field1_form_utils(W1, dt, deltac)
    else:
        raise ValueError(state)

    sol = Function(W1)
    problem = LinearVariationalProblem(a, L, sol, bcs=bcs, constant_jacobian=True)
    params = outer_parameters(patch, smooth_its, rtol, maxit, restart, smoother_ksp, mg_type, outer)
    _solver_counter[0] += 1
    prefix = f'inner{_solver_counter[0]}_'
    solver = LinearVariationalSolver(problem, solver_parameters=params, options_prefix=prefix)
    ksp = solver.snes.ksp
    ksp.setConvergenceHistory(length=maxit + 5, reset=True)

    t0 = time.perf_counter()
    solver.solve()
    t_first = time.perf_counter() - t0
    its = ksp.getIterationNumber()
    reason = ksp.getConvergedReason()
    hist = np.array(ksp.getConvergenceHistory()).copy()
    stats = collect_patch_stats(ksp)

    if view_file is not None and comm.rank == 0:
        viewer = PETSc.Viewer().createASCII(view_file, comm=comm)
        ksp.view(viewer)
        viewer.destroy()

    # ---- timing of the solve phase only (operators and factorisations are reused) ----
    t_solve = float('nan')
    for _ in range(timing_repeats):
        sol.assign(0)
        t0 = time.perf_counter()
        solver.solve()
        t1 = time.perf_counter() - t0
        t_solve = t1 if not (t_solve == t_solve) else min(t_solve, t1)
        assert ksp.getIterationNumber() == its, (ksp.getIterationNumber(), its)

    res0 = float(hist[0]) if len(hist) else float('nan')
    resn = float(hist[-1]) if len(hist) else float('nan')
    rate = (resn / res0) ** (1.0 / its) if its > 0 and res0 > 0 else float('nan')
    converged = (reason > 0) and (resn <= rtol * res0 * (1 + 1e-12))

    # ---- work proxies ----
    patch_solves = sum(s['napply'] * s['npatch'] for s in stats)
    work_patch = sum(s['napply'] * s['nz_patch_factors'] for s in stats)
    work_matmult = sum(s['napply'] * s['nz_level_mat'] for s in stats)   # one level MatMult per ASM application
    fine = stats[-1] if stats else {}
    row = dict(
        tag=tag, patch=patch, smooth_its=smooth_its, smoother_ksp=smoother_ksp, mg_type=mg_type, outer=outer,
        rhs=rhs, state=state, nranks=comm.size,
        nx=nx, ncol=ncol, nz=nz, Lx=Lx, Lz=Lz, dt=dt, K=K_eff, delta=delta, c=c, m=m, CM2=CM2,
        regime=regime(dt, Lx, Lz), dx=dxf, dz=dzf,
        c_over_dx2=c / dxf ** 2, c_over_dz2=c / dzf ** 2, c_over_m_dx2=c / (m * dxf ** 2),
        m_dz2_over_dx2=m * dzf ** 2 / dxf ** 2, c_pi2_over_m_Lz2=c * np.pi ** 2 / (m * Lz ** 2),
        ndofs=W1.dim(), its=its, converged=int(converged), reason=int(reason),
        res0=res0, res_final=resn, avg_rate=rate,
        patch_size=fine.get('patch_size_max', 0),
        npatch_levels='/'.join(str(s['npatch']) for s in stats),
        napply_levels='/'.join(str(s['napply']) for s in stats),
        asm_applies_per_mg='/'.join(f"{s['napply'] / max(its, 1):.2f}" for s in stats),
        patch_solves=patch_solves,
        nz_factor_per_patch=(fine.get('nz_patch_factors', 0) / max(fine.get('npatch', 1), 1)),
        work_patch_Mnz=work_patch / 1e6,
        work_matmult_Mnz=work_matmult / 1e6,
        work_total_Mnz=(work_patch + work_matmult) / 1e6,
        t_first=t_first, t_solve=t_solve,
    )
    solver = None
    return row, hist


# ----------------------------------------------------------------------------------
#  Command line
# ----------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description='Isolated field-1 multigrid solve (star vs vanka patch).',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tag', type=str, default='run')
    parser.add_argument('--out', type=str, default=None, help='CSV file (default results_<tag>.csv next to the script).')
    parser.add_argument('--nx', type=int, nargs='+', default=[20])
    parser.add_argument('--nz', type=int, nargs='+', default=[20])
    parser.add_argument('--Lx', type=float, nargs='+', default=[3.0e5])
    parser.add_argument('--Lz', type=float, nargs='+', default=[4.0e3])
    parser.add_argument('--dt', type=float, nargs='+', default=[50.0])
    parser.add_argument('--K', type=float, nargs='+', default=[1.0])
    parser.add_argument('--c', type=float, nargs='+', default=None, help='Prescribe c = dt/(2 delta) directly (overrides K).')
    parser.add_argument('--c_over_dx2', type=float, nargs='+', default=None, help='Prescribe c / dx^2 directly (overrides K and c).')
    parser.add_argument('--patch', type=str, nargs='+', default=['star', 'vanka'])
    parser.add_argument('--smooth_its', type=int, nargs='+', default=[6])
    parser.add_argument('--rhs', type=str, nargs='+', default=['random', 'bump'])
    parser.add_argument('--state', type=str, default='rest', choices=['rest', 'ic', 'rest_utils'])
    parser.add_argument('--refinement', type=int, default=2)
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--rtol', type=float, default=1.0e-8)
    parser.add_argument('--maxit', type=int, default=200)
    parser.add_argument('--restart', type=int, default=200)
    parser.add_argument('--smoother_ksp', type=str, default='gmres')
    parser.add_argument('--mg_type', type=str, default='full')
    parser.add_argument('--outer', type=str, default='fgmres')
    parser.add_argument('--timing_repeats', type=int, default=1)
    parser.add_argument('--view', action='store_true', help='Write the ksp_view of every solve to ksp_view_<tag>_<n>.txt')
    args, _ = parser.parse_known_args()

    here = os.path.dirname(os.path.abspath(__file__))
    out = args.out or os.path.join(here, f'results_{args.tag}.csv')
    hist_file = os.path.splitext(out)[0] + '_hist.jsonl'

    if args.c_over_dx2 is not None:
        shift_list = [('c_over_dx2', v) for v in args.c_over_dx2]
    elif args.c is not None:
        shift_list = [('c', v) for v in args.c]
    else:
        shift_list = [('K', v) for v in args.K]

    cache = {}
    n = 0
    # mesh parameters outermost so that a hierarchy is built once
    for nx, nz, Lx, Lz in itertools.product(args.nx, args.nz, args.Lx, args.Lz):
        for dt, (skind, sval), patch, sits, rhs in itertools.product(args.dt, shift_list, args.patch, args.smooth_its, args.rhs):
            kw = dict(K=None, c=None)
            if skind == 'K':
                kw['K'] = sval
            elif skind == 'c':
                kw['c'] = sval
            else:
                kw['c'] = sval * (Lx / (nx * 2 ** args.refinement)) ** 2
            n += 1
            view_file = os.path.join(here, f'ksp_view_{args.tag}_{n}_{patch}.txt') if args.view else None
            row, hist = run_case(cache, nx, nz, Lx, Lz, dt, patch=patch, smooth_its=sits, rhs=rhs,
                                 state=args.state, refinement=args.refinement, degree=args.degree,
                                 rtol=args.rtol, maxit=args.maxit, restart=args.restart,
                                 smoother_ksp=args.smoother_ksp, mg_type=args.mg_type, outer=args.outer,
                                 timing_repeats=args.timing_repeats, view_file=view_file, tag=args.tag, **kw)
            Print(f"[{args.tag} {n}] patch={patch:5s} sits={sits:2d} rhs={rhs:6s} nx={nx} nz={nz} Lz={Lz:g} dt={dt:g} "
                  f"K={row['K']:.4g} c/dx2={row['c_over_dx2']:.4g} m={row['m']:.4g} | its={row['its']:3d} "
                  f"conv={row['converged']} rate={row['avg_rate']:.3f} psolves={row['patch_solves']} "
                  f"work={row['work_total_Mnz']:.0f}Mnz t={row['t_solve']:.2f}s")
            if COMM_WORLD.rank == 0:
                new = not os.path.exists(out)
                with open(out, 'a', newline='') as f:
                    wr = csv.DictWriter(f, fieldnames=list(row.keys()))
                    if new:
                        wr.writeheader()
                    wr.writerow(row)
                with open(hist_file, 'a') as f:
                    rec = {k: row[k] for k in ('tag', 'patch', 'smooth_its', 'rhs', 'nx', 'nz', 'Lx', 'Lz', 'dt', 'K', 'c', 'its')}
                    rec['residuals'] = [float(v) for v in hist]
                    f.write(json.dumps(rec) + '\n')
            import gc
            gc.collect()


if __name__ == '__main__':
    main()
