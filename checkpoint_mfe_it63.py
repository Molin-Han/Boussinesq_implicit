from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils

distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 2)}
m = PeriodicIntervalMesh(10, 1.0,distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=2)
hierarchy = ExtrudedMeshHierarchy(mh, 1.0, layers=[10] * (3), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]

V = FunctionSpace(mesh, 'CG', 1)
x, y, z = SpatialCoordinate(mesh)
u = Function(V)
v = TestFunction(V)
f = Function(V).interpolate(x**2 + y**2 + z**2)
g = Constant(1.0)
equ = inner(grad(u), grad(v)) * dx - f * v * dx

params = {
    'snes_type': 'ksponly',
    'ksp_type':'gmres',
    'ksp_monitor':None,
    'ksp_view': ':mfe.txt',
    'pc_type': 'none',
}

nprob = NonlinearVariationalProblem(equ, u)
nsolver = NonlinearVariationalSolver(nprob, solver_parameters= params)

for i in range(128):
    if i == 0:
        with CheckpointFile('f.h5', 'w') as chk:
            print('Checkpointing the first solution.')
            chk.save_function(f, idx=i, name=f'sol_{i}')
            print('Checkpointed the first solution.')
    else:
        with CheckpointFile('f.h5', 'a') as chk:
            print(f'Checkpointing the solution at iteration {i}.')
            chk.save_function(Function(V).interpolate(f+g), idx=i, name=f'sol_{i}')
            print(f'Checkpointed the solution at iteration {i}.')

sol_it = Function(V, name='sol_it')

def monitor(ksp, iteration_number, norm0):
    print('The monitor starts to build the solution.')
    sol = ksp.buildSolution()
    print('The monitor finishes building the solution.')
    with sol_it.dat.vec_wo as it_vec:
        sol.copy(result=it_vec)
    if iteration_number == 0:
        with CheckpointFile('sol_its.h5', 'w') as chk:
            print('Checkpointing the first solution.')
            chk.save_function(sol_it, idx=iteration_number, name=f'sol_{iteration_number}')
            print('Checkpointed the first solution.')
    else:
        with CheckpointFile('sol_its.h5', 'a') as chk:
            print(f'Checkpointing the solution at iteration {iteration_number}.')
            chk.save_function(sol_it, idx=iteration_number, name=f'sol_{iteration_number}')
            print(f'Checkpointed the solution at iteration {iteration_number}.')

nsolver.snes.ksp.setMonitor(monitor)
nsolver.solve()