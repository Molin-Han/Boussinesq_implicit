from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


parser = ArgumentParser(
    description='MFE',
    formatter_class=ArgumentDefaultsHelpFormatter
)
parser.add_argument('--work', action='store_true', help='If True, the matrix is good, if False, a weird constant is added.')
args = parser.parse_known_args()
args = args[0]

nx = 3
nz = 3
length = 1.0
height = 1.0

m = PeriodicIntervalMesh(nx, length)
mesh = ExtrudedMesh(m, nz, layer_height=height/nz)

x, z = SpatialCoordinate(mesh)
deg = 1
CG = FiniteElement("CG", interval, deg)
DG = FiniteElement("DG", interval, deg-1)
CG_DG = TensorProductElement(CG, DG)
RT_horiz = HDivElement(CG_DG)
DG_CG = TensorProductElement(DG, CG)
RT_vert = HDivElement(DG_CG)
RT_e = RT_horiz + RT_vert
V = FunctionSpace(mesh, RT_e)

horiz = FiniteElement('DG', interval, deg-1)
vertical = FiniteElement('CG', interval, deg)
V_elt = TensorProductElement(horiz, vertical)
B = FunctionSpace(mesh, V_elt)
P = FunctionSpace(mesh, 'DG', deg-1)
f = Function(V).project(as_vector([sin(2 * pi * x), sin(2* pi * z / height)]))
g = Function(B).interpolate(sin(2*pi*x) + sin(2*pi*z))


if args.work:
    W = V * B
    U = Function(W)
    u, b = split(U)
    w, q = TestFunctions(W)
    F = inner(u, w) * dx + inner(div(u), div(w))*dx + b * q * dx - inner(f, w)*dx - g*q*dx
    params = {
        'mat_type': 'aij',
        'mat_view':':mat_noSchur_work_full.txt',
        'pc_type':'lu',
        'snes_monitor':None,
        'ksp_monitor':None,
        'pc_type':'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        'pc_fieldsplit_schur_precondition': 'full',
        'pc_fieldsplit_0_fields': '1',
        'pc_fieldsplit_1_fields': '0',
        'fieldsplit_0': {'pc_type':'lu','mat_view':':mat_noSchur_work_0.txt',},
        'fieldsplit_1': {'pc_type':'lu','mat_view':':mat_noSchur_work_1.txt',}
    }

else:
    W = V * B * P
    U = Function(W)
    u, b, p = split(U)
    w, q, phi = TestFunctions(W)

    F = inner(u, w) * dx + inner(div(u), div(w))*dx + b * q * dx - inner(f, w)*dx - g*q*dx
    F += p * phi * dx
    # F = inner(u, w) * dx - inner(p, div(w))*dx + b * q * dx - inner(f, w)*dx - g*q*dx
    # F += phi * div(u) * dx 
    params = {
        'mat_type': 'aij',
        'mat_view':':mat_noSchur_mixed_full.txt',
        'pc_type':'lu',
        'snes_monitor':None,
        'ksp_monitor':None,
        'pc_type':'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        'pc_fieldsplit_schur_precondition': 'full',
        'pc_fieldsplit_0_fields': '2',
        'pc_fieldsplit_1_fields': '0,1',
        'fieldsplit_0': {'pc_type':'lu','mat_view':':mat_noSchur_mixed_0.txt',},
        'fieldsplit_1': {'pc_type':'lu','mat_view':':mat_noSchur_mixed_1.txt',}
    }

nprob = NonlinearVariationalProblem(F, U)
nsolver = NonlinearVariationalSolver(nprob, solver_parameters=params)
nsolver.solve()