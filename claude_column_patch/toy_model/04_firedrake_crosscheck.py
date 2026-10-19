'''
Cross-check of the numpy toy model against Firedrake (serial, tiny meshes).

The same bilinear form a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w) is assembled by
Firedrake on a periodic extruded mesh with the element of utils.extrude_RT,

    HDiv(CG_k x DG_{k-1}) + HDiv(DG_{k-1} x CG_k),

and the additive Schwarz preconditioners are the ones of the solver script,

    star   : firedrake.ASMStarPC,  pc_star_construct_dim 0
    column : firedrake.ASMVankaPC, pc_vanka_construct_dim 1

with exact (LU) patch solves.  The matrix of the preconditioner B is obtained by applying the
PC to unit vectors, and the spectrum of B A is computed densely.

  degree 1 : the spectrum must coincide with the numpy toy model (rt0_model.py)
             -> this validates that the toy patches are the Firedrake patches.
  degree 2 : used as reference for the degree 2 toy model (rt1_model.py), if present.

usage: python 04_firedrake_crosscheck.py          (run inside the toy_model folder)
Output: 04_firedrake_crosscheck.txt, tables/firedrake_crosscheck_spectra.npz
'''
import os
import sys
import numpy as np
from firedrake import *
from firedrake.petsc import PETSc

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rt0_model as rt

HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, '04_firedrake_crosscheck.txt'), 'w')


def log(*args):
    s = ' '.join(str(a) for a in args)
    PETSc.Sys.Print(s)
    out.write(s + '\n')
    out.flush()


def firedrake_spectrum(nx, nz, Lx, Lz, mval, cval, degree, kind):
    '''Return (sorted eigenvalues of B A on the unconstrained DoFs, patch sizes, ndof).'''
    base = PeriodicIntervalMesh(nx, Lx)
    mesh = ExtrudedMesh(base, nz, layer_height=Lz / nz)
    CG = FiniteElement("CG", interval, degree)
    DG = FiniteElement("DG", interval, degree - 1)
    elt = HDivElement(TensorProductElement(CG, DG)) + HDivElement(TensorProductElement(DG, CG))
    V = FunctionSpace(mesh, elt)
    u = TrialFunction(V)
    w = TestFunction(V)
    mm = Constant(mval)
    cc = Constant(cval)
    a = (u[0] * w[0] + mm * u[1] * w[1] + cc * div(u) * div(w)) * dx
    bcs = [DirichletBC(V, as_vector([0., 0.]), "top"), DirichletBC(V, as_vector([0., 0.]), "bottom")]
    if kind == 'star':
        pc_params = {"pc_type": "python", "pc_python_type": "firedrake.ASMStarPC",
                     "pc_star_construct_dim": 0, "pc_star_sub_sub_pc_type": "lu"}
    else:
        pc_params = {"pc_type": "python", "pc_python_type": "firedrake.ASMVankaPC",
                     "pc_vanka_construct_dim": 1, "pc_vanka_sub_sub_pc_type": "lu"}
    params = {"mat_type": "aij", "ksp_type": "preonly", **pc_params}
    sol = Function(V)
    F = Cofunction(V.dual())
    problem = LinearVariationalProblem(a, F, sol, bcs=bcs)
    solver = LinearVariationalSolver(problem, solver_parameters=params)
    solver.solve()            # sets everything up (zero right-hand side)
    ksp = solver.snes.ksp
    pc = ksp.getPC()
    A, P = ksp.getOperators()
    n = A.getSize()[0]
    Ad = A.convert('dense').getDenseArray().copy()
    # matrix of the preconditioner
    x = A.createVecRight()
    y = A.createVecLeft()
    Bd = np.zeros((n, n))
    for j in range(n):
        x.set(0.0)
        x.setValue(j, 1.0)
        x.assemble()
        pc.apply(x, y)
        Bd[:, j] = y.getArray()
    # patch sizes
    asmpc = pc.getPythonContext().asmpc
    sizes = [ksp_.getOperators()[0].getSize()[0] for ksp_ in asmpc.getASMSubKSP()]
    # remove the constrained (boundary) DoFs: they are decoupled in A and in B
    bc_nodes = np.unique(np.concatenate([bc.nodes for bc in bcs]))
    free = np.setdiff1d(np.arange(n), bc_nodes)
    assert np.abs(Ad[np.ix_(free, bc_nodes)]).max() == 0.0
    assert np.abs(Bd[np.ix_(free, bc_nodes)]).max() < 1e-14 * np.abs(Bd).max()
    Af = Ad[np.ix_(free, free)]
    Bf = Bd[np.ix_(free, free)]
    assert np.allclose(Af, Af.T, rtol=1e-12, atol=1e-14 * np.abs(Af).max())
    sym_err = np.abs(Bf - Bf.T).max() / np.abs(Bf).max()
    L = np.linalg.cholesky(0.5 * (Af + Af.T))
    S = L.T @ (0.5 * (Bf + Bf.T)) @ L
    lam = np.sort(np.linalg.eigvalsh(0.5 * (S + S.T)))
    return lam, sizes, free.size, n, sym_err


Lx, Lz = 3.0e5, 4.0e3
results = {}
cases = [
    ('unit square cells, c = 1', 8, 6, 8.0, 6.0, 1.0, 1.0),
    ('flat cells, large c', 8, 6, 8.0, 1.8, 2.5, 1.0e4),
    ('nominal physical numbers (dt = 50 s, K = 1), 16 x 10', 16, 10, Lx, Lz,
     rt.physical_parameters(Lx, Lz, 50.0, 1.0)['m'], rt.physical_parameters(Lx, Lz, 50.0, 1.0)['c']),
    ('strongly stratified (dt = 1e5 s, K = 1), 16 x 10', 16, 10, Lx, Lz,
     rt.physical_parameters(Lx, Lz, 1e5, 1.0)['m'], rt.physical_parameters(Lx, Lz, 1e5, 1.0)['c']),
]

ok = True
log('degree 1: Firedrake (ASMStarPC / ASMVankaPC) against the numpy toy model')
for name, nx, nz, lx, lz, mval, cval in cases:
    for kind in ['star', 'column']:
        lam, sizes, nfree, n, sym_err = firedrake_spectrum(nx, nz, lx, lz, mval, cval, 1, kind)
        ltoy = rt.spectrum_2x2(nx, nz, lx / nx, lz / nz, mval, cval, kind)
        assert lam.size == ltoy.size, (lam.size, ltoy.size)
        err = np.max(np.abs(lam - ltoy) / ltoy)
        passed = err < 1e-6
        ok = ok and passed
        results[f'deg1|{name}|{kind}'] = lam
        log(f'  [{"PASS" if passed else "FAIL"}] {name:55s} {kind:6s} ndof {n} (free {nfree}), '
            f'Firedrake patch sizes {sorted(set(sizes))} (toy, without boundary DoFs: '
            f'{3 * nz - 2 if kind == "star" else 3 * nz - 1}), '
            f'lmin {lam[0]:.6e} lmax {lam[-1]:.6f} cond {lam[-1] / lam[0]:.4e}, '
            f'max rel. deviation from the toy {err:.1e}')

log('\ndegree 2: Firedrake spectra (reference for the degree 2 toy model)')
for name, nx, nz, lx, lz, mval, cval in cases:
    for kind in ['star', 'column']:
        lam, sizes, nfree, n, sym_err = firedrake_spectrum(nx, nz, lx, lz, mval, cval, 2, kind)
        results[f'deg2|{name}|{kind}'] = lam
        log(f'  {name:55s} {kind:6s} ndof {n} (free {nfree}), Firedrake patch sizes {sorted(set(sizes))}, '
            f'lmin {lam[0]:.6e} lmax {lam[-1]:.6f} cond {lam[-1] / lam[0]:.4e}, '
            f'eigenvalues below 0.1 lmax: {int(np.sum(lam < 0.1 * lam[-1]))}')

os.makedirs(os.path.join(HERE, 'tables'), exist_ok=True)
np.savez(os.path.join(HERE, 'tables', 'firedrake_crosscheck_spectra.npz'),
         case_names=np.array([c[0] for c in cases]),
         case_params=np.array([c[1:] for c in cases], dtype=float),
         **{k.replace('|', '__'): v for k, v in results.items()})
log('\nALL DEGREE 1 CROSS-CHECKS PASSED' if ok else '\nSOME CROSS-CHECKS FAILED')
out.close()
