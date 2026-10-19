#!/usr/bin/env python
'''
verify_against_solver.py
========================
Checks that inner_solve.py really isolates the multigrid of the solver script.

1. The UNMODIFIED solver script Nonlinear_Boussinesq_Irk_SC.py is executed in this
   process (exec of its source, one time step, the SNES is stopped after the first
   Newton step through -firedrake_1_snes_rtol so that the Jacobian is the one
   linearised about the initial condition: rest + buoyancy perturbation, stage = 0).
   It runs inside scratch/solver_run_<patch>/ because it writes .pvd/.vtu files.
2. The fieldsplit_1 operators (S, A11) and the multigrid PC object of the solver are
   taken out of the PETSc solver hierarchy.
3. Checks:
   (a) S x = A11 x for random x (the Schur complement of the auxiliary operator IS A11,
       because the pressure row of the auxiliary operator does not depend on u);
   (b) A11 of the solver == matrix assembled by inner_solve.field1_form_utils (state 'ic');
   (c) the hand-written form linearised about rest (inner_solve.field1_form_rest) ==
       field1_form_utils at rest, and how far the matrix at 'ic' is from the one at rest;
   (d) FGMRES(rtol 1e-8) preconditioned by the solver's OWN multigrid object gives the
       same iteration counts / residual histories as inner_solve.run_case.

usage:  python verify_against_solver.py --patch vanka [--K 1 --dt 50 --nx 20 --nz 20]
'''
import os
import sys
import json

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('BOUSSINESQ_REPO', '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
SOLVER = os.path.join(REPO, 'Nonlinear_Boussinesq_Irk_SC.py')


def _own_args(argv):
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('--patch', default='vanka', choices=['star', 'vanka'])
    p.add_argument('--K', type=float, default=1.0)
    p.add_argument('--dt', type=float, default=50.0)
    p.add_argument('--nx', type=int, default=20)
    p.add_argument('--nz', type=int, default=20)
    p.add_argument('--Lx', type=float, default=3.0e5)
    p.add_argument('--Lz', type=float, default=4.0e3)
    a, _ = p.parse_known_args(argv)
    return a


_own = _own_args(sys.argv[1:])
_m = 1.0 + (_own.dt * 1.0e-2) ** 2 / 4.0
import math                                                        # noqa: E402
_CM2 = max((_own.Lx / (2 * math.pi)) ** 2, _m * (_own.Lz / math.pi) ** 2)
_shift = _own.K * _own.dt / (2.0 * _CM2)

_workdir = os.path.join(HERE, 'scratch', f'solver_run_{_own.patch}')
os.makedirs(_workdir, exist_ok=True)
os.chdir(_workdir)
for _p in (REPO, HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)

# command line seen by the solver script (argparse) and by PETSc
sys.argv = [SOLVER,
            '--nx', str(_own.nx), '--nz', str(_own.nz),
            '--length', repr(_own.Lx), '--height', repr(_own.Lz),
            '--dt', repr(_own.dt), '--tmax', repr(_own.dt),
            '--shift', repr(_shift), '--maxit', '100', '--no_rotation',
            '--patch', _own.patch,
            '-firedrake_1_snes_rtol', '1e-3']

with open(SOLVER) as _f:
    _src = _f.read()
exec(compile(_src, SOLVER, 'exec'), globals())      # runs the solver script, defines `stepper`

# ------------------------------------------------------------------------------------
#  From here on: our own code.  (The exec above defined many names, among them
#  print = PETSc.Sys.Print, stepper, W, args.)
# ------------------------------------------------------------------------------------
import numpy as _np                                                # noqa: E402
import firedrake as _fd                                            # noqa: E402
import firedrake.dmhooks as _dmhooks                               # noqa: E402
from firedrake.petsc import PETSc as _PETSc                        # noqa: E402
import inner_solve as _ins                                         # noqa: E402
import claude_paper_shift as _cps                                  # noqa: E402

_report = {'patch': _own.patch, 'K': _own.K, 'dt': _own.dt, 'nx': _own.nx, 'nz': _own.nz,
           'Lx': _own.Lx, 'Lz': _own.Lz, 'shift': _shift}

# shift formula of inner_solve.py == claude_paper_shift.py
_report['shift_formula_relative_difference'] = abs(
    _ins.shift_threshold(_own.dt, _own.Lx, _own.Lz) - _cps.shift_threshold(_own.dt, _own.Lx, _own.Lz)
) / _cps.shift_threshold(_own.dt, _own.Lx, _own.Lz)

_snes = stepper.solver.snes                                         # noqa: F821
_report['solver_options_prefix'] = stepper.solver.options_prefix   # noqa: F821
_report['solver_newton_iterations'] = int(_snes.getIterationNumber())
_report['solver_outer_fgmres_iterations_first_newton_step'] = int(_snes.getLinearSolveIterations())
assert _snes.getIterationNumber() == 1, 'the Jacobian must be the one of the first Newton step'

_pc_py = _snes.ksp.getPC().getPythonContext()                       # HDivSchurPC
_fs = _pc_py.pc                                                     # PCFIELDSPLIT
assert _fs.getType() == 'fieldsplit'
_ksp0, _ksp1 = _fs.getFieldSplitSubKSP()
_S, _A11 = _ksp1.getOperators()
_mg = _ksp1.getPC()
assert _mg.getType() == 'mg'
_report['solver_A11_size'] = int(_A11.getSize()[0])
_report['solver_A11_nnz'] = float(_A11.getInfo()['nz_used'])
_report['solver_fieldsplit1_Amat_type'] = _S.getType()

# ---- (a) S == A11 --------------------------------------------------------------------
_x = _A11.createVecRight()
_x.setRandom()
_y1 = _A11.createVecLeft()
_y2 = _A11.createVecLeft()
_A11.mult(_x, _y1)
_S.mult(_x, _y2)
_y2.axpy(-1.0, _y1)
_report['rel_diff_Schur_vs_A11_on_random_vector'] = _y2.norm() / _y1.norm()

# ---- standalone operators on a second, identical hierarchy -------------------------------
_mh = _ins.build_hierarchy(_own.nx, _own.nz, _own.Lx, _own.Lz, 2)
_W1 = _ins.build_space(_mh[-1], 2)
_delta = _fd.Constant(_shift)
_state_ic = _fd.Function(_W1)
_state_ic.subfunctions[2].assign(_ins.initial_buoyancy(_W1, _own.Lx, _own.Lz))
_J_ic, _bcs_ic, _, _ = _ins.field1_form_utils(_W1, _own.dt, _delta, state=_state_ic)
_J_ru, _bcs_ru, _, _ = _ins.field1_form_utils(_W1, _own.dt, _delta)
_a_rest = _ins.field1_form_rest(_W1, _fd.Constant(_own.dt), _delta)
_bcs_rest = _ins.boundary_conditions(_W1)
_M_ic = _fd.assemble(_J_ic, bcs=_bcs_ic, mat_type='aij').petscmat
_M_ru = _fd.assemble(_J_ru, bcs=_bcs_ru, mat_type='aij').petscmat
_M_rest = _fd.assemble(_a_rest, bcs=_bcs_rest, mat_type='aij').petscmat


def _rel_diff(A, B):
    D = A.copy()
    D.axpy(-1.0, B, structure=_PETSc.Mat.Structure.DIFFERENT_NONZERO_PATTERN)
    return D.norm(_PETSc.NormType.FROBENIUS) / B.norm(_PETSc.NormType.FROBENIUS)


# ---- (b), (c) ----------------------------------------------------------------------------
_report['rel_diff_solverA11_vs_standalone_ic'] = _rel_diff(_A11, _M_ic)
_report['rel_diff_handwritten_rest_vs_utils_rest'] = _rel_diff(_M_rest, _M_ru)
_report['rel_diff_ic_vs_rest'] = _rel_diff(_M_ic, _M_ru)
_report['standalone_nnz_utils_form'] = float(_M_ic.getInfo()['nz_used'])
_report['standalone_nnz_handwritten_form'] = float(_M_rest.getInfo()['nz_used'])

# ---- (d) the solver's own multigrid object inside our FGMRES -----------------------------
_cache = {}
_results = []
for _rhs_kind in ('random', 'bump'):
    # standalone (inner_solve.py), both states
    _rows = {}
    for _state in ('ic', 'rest_utils', 'rest'):
        _row, _hist = _ins.run_case(_cache, _own.nx, _own.nz, _own.Lx, _own.Lz, _own.dt, K=_own.K,
                                    patch=_own.patch, rhs=_rhs_kind, state=_state, timing_repeats=0)
        _rows[_state] = (_row['its'], [float(v) for v in _hist])
    # right-hand side vector with the boundary rows zeroed, in the ordering of A11
    _key = list(_cache.keys())[0]
    _L = _cache[_key]['rhs'][_rhs_kind]
    _Wc = _cache[_key]['W1']
    _Lc = _L.copy(deepcopy=True)
    for _bc in _ins.boundary_conditions(_Wc):
        _Lc.dat[0].data[_bc.nodes] = 0.0
    _b = _A11.createVecLeft()
    with _Lc.dat.vec_ro as _v:
        _b.setArray(_v.getArray().copy())
    _xs = _A11.createVecRight()
    _ksp = _PETSc.KSP().create(comm=_A11.comm)
    _ksp.setOptionsPrefix('verify_outer_')
    _ksp.setType('fgmres')
    _ksp.setPCSide(_PETSc.PC.Side.RIGHT)
    _ksp.setTolerances(rtol=1.0e-8, atol=0.0, max_it=200)
    _ksp.setGMRESRestart(200)
    _ksp.setPC(_mg)                       # the multigrid object of the solver script, operators (S, A11)
    _ksp.setConvergenceHistory(length=210, reset=True)
    _dm = _ksp1.getDM()
    with _dmhooks.add_hooks(_snes.ksp.getPC().getDM(), _pc_py, appctx=_pc_py._ctx_ref):
        _ksp.solve(_b, _xs)
    _h = [float(v) for v in _ksp.getConvergenceHistory()]
    _res = dict(rhs=_rhs_kind,
                its_solver_mg_object=int(_ksp.getIterationNumber()),
                its_standalone_ic=_rows['ic'][0],
                its_standalone_rest_utils=_rows['rest_utils'][0],
                its_standalone_rest_handwritten=_rows['rest'][0],
                max_rel_diff_history_solver_vs_standalone_ic=float(
                    max(abs(a - b) / max(abs(b), 1e-300) for a, b in zip(_h, _rows['ic'][1]))
                ) if len(_h) == len(_rows['ic'][1]) else None,
                history_solver_mg_object=_h,
                history_standalone_ic=_rows['ic'][1],
                history_standalone_rest=_rows['rest'][1])
    _results.append(_res)
    _ksp.destroy()
_report['inner_solves'] = _results

_out = os.path.join(HERE, f'verification_{_own.patch}.json')
with open(_out, 'w') as _f:
    json.dump(_report, _f, indent=1)
for _k, _v in _report.items():
    if _k != 'inner_solves':
        _PETSc.Sys.Print(f'{_k}: {_v}')
for _r in _results:
    _PETSc.Sys.Print({k: v for k, v in _r.items() if not k.startswith('history')})
