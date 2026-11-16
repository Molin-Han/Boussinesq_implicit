from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print

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

W = V * B * P
U = Function(W)
u, b, p = split(U)
w, q, phi = TestFunctions(W)

f = Function(V).project(as_vector([sin(2 * pi * x), sin(2* pi * z / height)]))
g = Function(B).interpolate(sin(2*pi*x) + sin(2*pi*z))

F = inner(u, w) * dx + inner(div(u), div(w))*dx + b * q * dx - inner(f, w)*dx - g*q*dx
# F += phi * div(u) * dx 
F += p * phi * dx

params_mixed = {
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

params_block = {
    'mat_type': 'aij',
    'mat_view':':mat_noSchur_block_full.txt',
    'pc_type':'lu',
    'snes_monitor':None,
    'ksp_monitor':None,
    'pc_type':'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_fact_type': 'full',
    'pc_fieldsplit_schur_precondition': 'full',
    'pc_fieldsplit_0_fields': '2',
    'pc_fieldsplit_1_fields': '0,1',
    'fieldsplit_0': {'pc_type':'lu','mat_view':':mat_noSchur_block_0.txt',},
    'fieldsplit_1': {'pc_type':'fieldsplit',
        'pc_fieldsplit_type':'additive',
        'fieldsplit':{'pc_type':'lu'},
        'fieldsplit_0_mat_view':':mat_noSchur_block_1.txt',
        'fieldsplit_1_mat_view':':mat_noSchur_block_2.txt',
        }
}


nprob = NonlinearVariationalProblem(F, U)
# nsolver = NonlinearVariationalSolver(nprob, solver_parameters=params_mixed)
nsolver = NonlinearVariationalSolver(nprob, solver_parameters=params_block)
nsolver.solve()