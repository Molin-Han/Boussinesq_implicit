"""
Item 3/4, operator level: is the preconditioner of a parallel run the same operator as that of
a serial run (and of a second run with the same number of ranks)?

Iteration counts are a blunt instrument for this question (the smoother is GMRES, i.e. nonlinear,
and the convergence test is a threshold), so the pieces of the preconditioner are applied to
fixed vectors and the results are stored DoF by DoF under partition-independent keys:

    P x           assembled (u, u_y, b) block of the auxiliary operator on every level
    ASM x         the additive Schwarz operator alone (linear) on every level with smoother
    smoother x    6 GMRES iterations preconditioned by ASM, zero initial guess (nonlinear in x)
    coarse x      LU on the coarsest level
    MG x          one application of the complete fieldsplit_1 preconditioner (full MG)

How.  The real solver script is executed for one time step (exec in __main__).  Then the state is
reset to a partition-independent one (velocity 0, b = nodal interpolant of the initial profile,
stages 0) and the time step is taken again with a KSP monitor that does the work at iteration 0
of the first linear solve, i.e. with the preconditioner set up at that state and before it has
been used.  The test vectors are functions of the keys only.

Usage (each run writes runs/operator_<tag>/outputs.npz and report.txt):
    python 03_operator_serial_vs_parallel.py --patch vanka --nx 4 --nz 3 --tag vanka_small_np1
    mpiexec --bind-to none -n 2 python 03_operator_serial_vs_parallel.py --patch vanka --nx 4 --nz 3 --tag vanka_small_np2
    python 03b_compare_operator_outputs.py runs/operator_vanka_small_np1 runs/operator_vanka_small_np2
"""
import os
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--patch', choices=['star', 'vanka'], default='vanka')
parser.add_argument('--nx', type=int, default=4)
parser.add_argument('--nz', type=int, default=3)
parser.add_argument('--degree', type=int, default=2)
parser.add_argument('--refinement', type=int, default=2)
parser.add_argument('--length', type=float, default=300000.0)
parser.add_argument('--height', type=float, default=4000.0)
parser.add_argument('--dt', type=float, default=50.0)
parser.add_argument('--shift', type=float, default=1.096622711232e-08)
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

_comm = MPI.COMM_WORLD
_rank, _size = _comm.rank, _comm.size
_L = float(my_args.length)
RUN = os.path.join(HERE, "runs", "operator_" + my_args.tag)
if _rank == 0:
    os.makedirs(RUN, exist_ok=True)
_comm.Barrier()
os.chdir(RUN)

lines = []
results = {}          # name -> (keys..., values) gathered on rank 0
state = {'done': False}


def out(s=""):
    lines.append(s)


def gather_sorted(T, vec):
    """Owned entries of vec with their keys, gathered on rank 0 and sorted by key."""
    V_idx = np.nonzero(T['owned'])[0]
    lg = T['lgmap_indices']
    rstart, rend = vec.getOwnershipRange()
    g = lg[V_idx]
    assert ((g >= rstart) & (g < rend)).all()
    vals = vec.getArray(readonly=True)[g - rstart].copy()
    f, c, bx, pos = pt.key_arrays(T, V_idx)
    allv = _comm.gather((f, c, bx, pos, vals), root=0)
    if _rank != 0:
        return None
    f = np.concatenate([a[0] for a in allv])
    c = np.concatenate([a[1] for a in allv])
    bx = np.concatenate([a[2] for a in allv])
    pos = np.concatenate([a[3] for a in allv])
    vals = np.concatenate([a[4] for a in allv])
    order = np.lexsort((pos, bx, c, f))
    return f[order], c[order], bx[order], pos[order], vals[order]


def fill(T, vec, which):
    V_idx = np.nonzero(T['owned'])[0]
    lg = T['lgmap_indices']
    rstart, rend = vec.getOwnershipRange()
    arr = vec.getArray()
    arr[lg[V_idx] - rstart] = pt.key_values(T, V_idx, which)
    vec.assemble()


def store(name, T, vec):
    res = gather_sorted(T, vec)
    nrm = vec.norm()
    if _rank == 0:
        results[name] = res
        out(f"   {name:34s} ||y||_2 = {nrm:.15e}")


def analyse(ksp):
    pc = ksp.getPC()
    ctx = pc.getPythonContext()
    inner = ctx.pc
    assert inner.getType() == 'fieldsplit', inner.getType()
    sub = inner.getFieldSplitSubKSP()
    assert len(sub) == 2
    ksp1 = sub[1]
    mg = ksp1.getPC()
    assert mg.getType() == 'mg', mg.getType()
    nlev = mg.getMGLevels()
    out(f"MG levels: {nlev}")
    for lev in range(nlev):
        if lev == 0:
            k = mg.getMGCoarseSolve()
        else:
            k = mg.getMGSmoother(lev)
        V = get_function_space(k.getDM())
        T = pt.describe_space(V, _L)
        T['lgmap_indices'] = V.dof_dset.lgmap.indices.copy()
        A, P = k.getOperators()
        n_glob = P.getSize()[0]
        dm = T['mesh'].topology_dm
        cs, ce = dm.getHeightStratum(0)
        ps, _pe = dm.getChart()
        owned_cells = sorted(pt.canon(T['px'][c - ps], _L) for c in range(cs, ce) if not T['ghost'][c - ps])
        allc = _comm.gather(owned_cells, root=0)
        if _rank == 0:
            ncols = sum(len(a) for a in allc)
            dx = _L / ncols
            out("")
            out(f"level {lev}: {ncols} base cells, {n_glob} DoFs, operator type {P.getType()}, "
                f"solver {k.getType()} + {k.getPC().getType()}")
            for r, a in enumerate(allc):
                idxs = sorted(int(round((xc - 0.5 * dx) / dx)) for xc in a)
                runs_ = []
                for i in idxs:
                    if runs_ and i == runs_[-1][1] + 1:
                        runs_[-1][1] = i
                    else:
                        runs_.append([i, i])
                out(f"   rank {r} owns the base cells number {', '.join(f'{a_}..{b_}' for a_, b_ in runs_)}")
            results[f"partition_level{lev}"] = (np.array([len(a) for a in allc]),)
        for which in ('smooth', 'rough'):
            x = P.createVecRight()
            y = P.createVecLeft()
            fill(T, x, which)
            store(f"L{lev} x[{which}]", T, x)
            P.mult(x, y)
            store(f"L{lev} P x[{which}]", T, y)
            if lev == 0:
                y.set(0.0)
                k.solve(x, y)
                store(f"L{lev} coarse x[{which}]", T, y)
            else:
                spc = k.getPC()
                actx = spc.getPythonContext()
                asmpc = actx.asmpc
                y.set(0.0)
                asmpc.apply(x, y)
                store(f"L{lev} ASM x[{which}]", T, y)
                y.set(0.0)
                k.solve(x, y)
                store(f"L{lev} smoother x[{which}]", T, y)
                # the same smoother a second time on the same input: is one process reproducible?
                y2 = y.duplicate()
                y2.set(0.0)
                k.solve(x, y2)
                y2.axpy(-1.0, y)
                d = y2.norm() / y.norm()
                if _rank == 0:
                    out(f"   {'   (second application, same run)':34s} rel. difference {d:.2e}")
                if lev == nlev - 1:
                    y.set(0.0)
                    ksp1.solve(x, y)
                    store(f"L{lev} MG x[{which}]", T, y)
                    y2.set(0.0)
                    ksp1.solve(x, y2)
                    y2.axpy(-1.0, y)
                    d = y2.norm() / y.norm()
                    if _rank == 0:
                        out(f"   {'   (second application, same run)':34s} rel. difference {d:.2e}")
            x.destroy()
            y.destroy()


def my_monitor(ksp, it, rnorm):
    if not state['done']:
        state['done'] = True
        analyse(ksp)


sys.argv = [SOLVER,
            '--nx', str(my_args.nx), '--nz', str(my_args.nz), '--degree', str(my_args.degree),
            '--refinement', str(my_args.refinement),
            '--length', repr(my_args.length), '--height', repr(my_args.height),
            '--dt', repr(my_args.dt), '--tmax', repr(my_args.dt),          # one time step
            '--shift', repr(my_args.shift), '--maxit', '100', '--no_rotation', '--patch', my_args.patch,
            '--timing']                                                    # no VTK output, no monitors
src = open(SOLVER).read()
_saved = dict(globals())
exec(compile(src, SOLVER, 'exec'), globals())
_solver = dict(globals())
globals().pop('print', None)
globals().update(_saved)

# ---- partition-independent state, then the same time step again with the monitor ----------------
stepper = _solver['stepper']
U = _solver['U']
t = _solver['t']
mesh = _solver['mesh']
for s in U.subfunctions:
    s.assign(0.0)
X = fd.SpatialCoordinate(mesh)
b_sub = U.subfunctions[2]
b_sub.interpolate(3e-4 * fd.sin(fd.pi * X[2] / my_args.height) / (1 + ((X[0] - my_args.length / 2) ** 2) / 5000.0 ** 2))
stepper.stages.assign(0.0)
t.assign(0.0)
out(f"patch = {my_args.patch}, nx = {my_args.nx}, nz = {my_args.nz}, degree = {my_args.degree}, "
    f"dt = {my_args.dt}, shift = {my_args.shift}, MPI ranks = {_size}")
out(f"state for the comparison: u = 0, u_y = 0, p = 0, b = nodal interpolant, ||b||_l2 = "
    f"{np.sqrt(_comm.allreduce(float((b_sub.dat.data_ro ** 2).sum()), op=MPI.SUM)):.15e}")
stepper.solver.snes.ksp.setMonitor(my_monitor)
stepper.advance()
its = stepper.solver.snes.getLinearSolveIterations()
out("")
out(f"after the analysis the time step was completed: Newton iterations "
    f"{stepper.solver.snes.getIterationNumber()}, cumulative FGMRES iterations {its}")

if _rank == 0:
    save = {}
    for name, res in results.items():
        key = name.replace(' ', '_').replace('[', '_').replace(']', '')
        if name.startswith('partition'):
            save[key] = res[0]
        else:
            f, c, bx, pos, vals = res
            save[key + "__vals"] = vals
            save[key + "__keys"] = np.vstack([f, c, np.round(bx * 1000).astype(np.int64), pos])
    np.savez_compressed("outputs.npz", **save)
    text = "\n".join(lines)
    sys.stdout.write(text + "\n")
    sys.stdout.flush()
    with open("report.txt", "w") as fh:
        fh.write(text + "\n")
