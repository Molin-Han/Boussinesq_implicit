"""
Item 5 (fairness): how good is ONE multigrid cycle with star / vanka smoother as a solver for the
(u, u_y, b) block, independently of the right-hand side of the test problem?

In the solver script the quality of the smoother is only seen through the outer FGMRES count, for
one particular (smooth) right-hand side.  Here the inner problem

        S y = x ,      S = A11 - A10 inv(A00) A01  (= A11, the block with (1/delta) div-div)

is solved by FGMRES preconditioned by exactly the multigrid preconditioner object of the solver
script (pc_mg_type full, GMRES(6)+ASM smoother, LU on the coarsest level), for
    x = smooth   (smooth function of x and z)
    x = rough    (pseudo-random, excites every mode)
and the number of iterations for a residual reduction of 1e-6 and 1e-10 is recorded, together
with the mean convergence factor.  The objects are taken from the real solver at iteration 0 of
the first linear solve of a time step that starts from a partition-independent state.

Several (dt, K) can be done in one process: dt and the shift are Constants, they are re-assigned
and the time step is started again (it is stopped after the first Newton step).

Usage:
    mpiexec --bind-to none -n 2 python 05_inner_solver_test.py --patch vanka --nx 20 --nz 20 \
          --dts 50 5000 --Ks 0.1 1 10 --tag vanka_nx20
Output: runs/inner_<tag>/report.txt, results.json
"""
import os
import sys
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--patch', choices=['star', 'vanka'], default='vanka')
parser.add_argument('--nx', type=int, default=20)
parser.add_argument('--nz', type=int, default=20)
parser.add_argument('--degree', type=int, default=2)
parser.add_argument('--refinement', type=int, default=2)
parser.add_argument('--length', type=float, default=300000.0)
parser.add_argument('--height', type=float, default=4000.0)
parser.add_argument('--dts', type=float, nargs='+', default=[50.0])
parser.add_argument('--Ks', type=float, nargs='+', default=[1.0])
parser.add_argument('--maxit', type=int, default=150)
parser.add_argument('--tag', type=str, required=True)
my_args, _ = parser.parse_known_args()

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/home/mh2520/molin/phd/yr2/Boussinesq_implicit"
SOLVER = os.path.join(REPO, "Nonlinear_Boussinesq_Irk_SC.py")
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import numpy as np                                           # noqa: E402
from mpi4py import MPI                                       # noqa: E402
import firedrake as fd                                       # noqa: E402
from firedrake.dmhooks import get_function_space             # noqa: E402
from firedrake.petsc import PETSc                            # noqa: E402
import patch_tools as pt                                     # noqa: E402
import claude_paper_shift as cps                             # noqa: E402  (the shift formula of the sweeps)

_comm = MPI.COMM_WORLD
_rank, _size = _comm.rank, _comm.size
_L = float(my_args.length)
RUN = os.path.join(HERE, "runs", "inner_" + my_args.tag)
if _rank == 0:
    os.makedirs(RUN, exist_ok=True)
_comm.Barrier()
os.chdir(RUN)

lines = []
results = []
state = {'armed': False, 'combo': None}


def out(s=""):
    lines.append(s)
    if _rank == 0:
        sys.stdout.write(s + "\n")
        sys.stdout.flush()


def fill(T, vec, which):
    V_idx = np.nonzero(T['owned'])[0]
    lg = T['lgmap_indices']
    rstart, rend = vec.getOwnershipRange()
    arr = vec.getArray()
    arr[lg[V_idx] - rstart] = pt.key_values(T, V_idx, which)
    vec.assemble()


_tables = {}


def analyse(ksp):
    pc = ksp.getPC()
    inner = pc.getPythonContext().pc
    sub = inner.getFieldSplitSubKSP()
    ksp1 = sub[1]
    mg = ksp1.getPC()
    assert mg.getType() == 'mg'
    nlev = mg.getMGLevels()
    S, P = ksp1.getOperators()
    if 'T' not in _tables:
        V = get_function_space(ksp1.getDM())
        T = pt.describe_space(V, _L)
        T['lgmap_indices'] = V.dof_dset.lgmap.indices.copy()
        _tables['T'] = T
    T = _tables['T']
    k = PETSc.KSP().create(comm=ksp1.comm)
    k.setOperators(S, P)
    k.setType('fgmres')
    k.setGMRESRestart(my_args.maxit)
    k.setPCSide(PETSc.PC.Side.RIGHT)
    k.setPC(mg)
    k.setTolerances(rtol=1e-10, atol=1e-300, divtol=1e10, max_it=my_args.maxit)
    k.setConvergenceHistory()
    k.setInitialGuessNonzero(False)
    rec = dict(state['combo'])
    rec.update(patch=my_args.patch, nx=my_args.nx, nz=my_args.nz, height=my_args.height, nranks=_size,
               ndof=int(P.getSize()[0]), levels=int(nlev))
    for which in ('smooth', 'rough'):
        x = P.createVecRight()
        y = P.createVecLeft()
        fill(T, x, which)
        y.set(0.0)
        k.solve(x, y)
        h = np.array(k.getConvergenceHistory())
        rel = h / h[0]
        it6 = int(np.argmax(rel <= 1e-6)) if (rel <= 1e-6).any() else -1
        it10 = int(np.argmax(rel <= 1e-10)) if (rel <= 1e-10).any() else -1
        n = len(rel) - 1
        rate_all = float(rel[-1] ** (1.0 / n)) if n > 0 else float('nan')
        # asymptotic factor: geometric mean over the second half of the iterations
        half = max(1, n // 2)
        rate_late = float((rel[-1] / rel[half]) ** (1.0 / (n - half))) if n > half else float('nan')
        rec[which] = dict(it_1e6=it6, it_1e10=it10, its=n, reason=int(k.getConvergedReason()),
                          rate_mean=rate_all, rate_second_half=rate_late,
                          history=[float(v) for v in rel])
        x.destroy()
        y.destroy()
    k.setPC(PETSc.PC().create(comm=ksp1.comm))      # detach the multigrid PC before the KSP is destroyed
    k.destroy()
    results.append(rec)
    out(f"   dt = {rec['dt']:9g}  K = {rec['K']:7g}  shift = {rec['shift']:.6e}  c/dx^2 = {rec['c_over_dx2']:9.3g}  "
        f"c/dz^2 = {rec['c_over_dz2']:9.3g}  m = {rec['m']:8.4g} | smooth rhs: its(1e-6) {rec['smooth']['it_1e6']:3d}, "
        f"its(1e-10) {rec['smooth']['it_1e10']:3d}, factor {rec['smooth']['rate_second_half']:.3f} | "
        f"rough rhs: its(1e-6) {rec['rough']['it_1e6']:3d}, its(1e-10) {rec['rough']['it_1e10']:3d}, "
        f"factor {rec['rough']['rate_second_half']:.3f}")


def my_monitor(ksp, it, rnorm):
    if state['armed']:
        state['armed'] = False
        analyse(ksp)
        # stop the time step as early as possible: one FGMRES iteration, one Newton step
        ksp.setTolerances(max_it=1)
        _solver['stepper'].solver.snes.setTolerances(max_it=1)


dt0 = my_args.dts[0]
sys.argv = [SOLVER,
            '--nx', str(my_args.nx), '--nz', str(my_args.nz), '--degree', str(my_args.degree),
            '--refinement', str(my_args.refinement),
            '--length', repr(my_args.length), '--height', repr(my_args.height),
            '--dt', repr(dt0), '--tmax', repr(dt0),
            '--shift', repr(cps.shift_value(1.0, dt0, my_args.length, my_args.height)),
            '--maxit', '100', '--no_rotation', '--patch', my_args.patch, '--timing']
src = open(SOLVER).read()
_saved = dict(globals())
exec(compile(src, SOLVER, 'exec'), globals())
_solver = dict(globals())
globals().pop('print', None)
globals().update(_saved)

stepper = _solver['stepper']
U = _solver['U']
t = _solver['t']
mesh = _solver['mesh']
dtc = _solver['dt']
shiftc = _solver['shift']
X = fd.SpatialCoordinate(mesh)
ncols = my_args.nx * 2 ** my_args.refinement
dx = my_args.length / ncols
dz = my_args.height / my_args.nz
out(f"patch = {my_args.patch}, columns on the finest level = {ncols} (dx = {dx:g} m), layers = {my_args.nz} "
    f"(dz = {dz:g} m), dx/dz = {dx / dz:g}, L_z/L_x = {my_args.height / my_args.length:g}, MPI ranks = {_size}")
stepper.solver.snes.ksp.setMonitor(my_monitor)
for dtv in my_args.dts:
    for K in my_args.Ks:
        sh = cps.shift_value(K, dtv, my_args.length, my_args.height)
        c = dtv / (2.0 * sh)
        state['combo'] = dict(dt=dtv, K=K, shift=sh, c=c, c_over_dx2=c / dx ** 2, c_over_dz2=c / dz ** 2,
                              m=cps.m_factor(dtv), regime=cps.regime(dtv, my_args.length, my_args.height))
        dtc.assign(dtv)
        shiftc.assign(sh)
        for s in U.subfunctions:
            s.assign(0.0)
        U.subfunctions[2].interpolate(3e-4 * fd.sin(fd.pi * X[2] / my_args.height)
                                      / (1 + ((X[0] - my_args.length / 2) ** 2) / 5000.0 ** 2))
        stepper.stages.assign(0.0)
        t.assign(0.0)
        stepper.solver.snes.ksp.setTolerances(max_it=100)
        stepper.solver.snes.setTolerances(max_it=10)
        state['armed'] = True
        try:
            stepper.advance()
        except fd.ConvergenceError:
            pass
        assert not state['armed'], "the monitor was not called"

if _rank == 0:
    with open("results.json", "w") as fh:
        json.dump(results, fh)
    with open("report.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")
