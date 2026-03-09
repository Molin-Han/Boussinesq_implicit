from firedrake import *
from petsc4py import PETSc
print = PETSc.Sys.Print

# m = PeriodicIntervalMesh(10, 1.0)
m = UnitIntervalMesh(10)
mesh = ExtrudedMesh(m, 20, extrusion_type='uniform')
x, y = SpatialCoordinate(mesh)
CG = FiniteElement("CG", interval, 1)
DG = FiniteElement("DG", interval, 0)
CG_DG = TensorProductElement(CG, DG)
RT_horiz = HDivElement(CG_DG)
DG_CG = TensorProductElement(DG, CG)
RT_vert = HDivElement(DG_CG)
RT_e = RT_horiz + RT_vert
V_2D = FunctionSpace(mesh, RT_e)
Pressure = FunctionSpace(mesh, 'DG', 0)
W = V_2D * Pressure
U = Function(W)
u, p = split(U)
w, phi = TestFunctions(W)

f = Function(V_2D).project(as_vector([x**2 + y**2, x*y]))
equ = inner(u, w) * dx  + p * phi * dx - inner(f, w) * dx

params = {
    'snes_type': 'ksponly',
    'snes_max_it':200,
    'ksp_type':'gmres',
    'ksp_monitor':None,
    'ksp_atol': 0,
    'ksp_rtol':1e-50,
    'ksp_max_it': 200,
    'ksp_gmres_restart':1000,
    'pc_type': 'none',
}

nprob = NonlinearVariationalProblem(equ, U)
nsolver = NonlinearVariationalSolver(nprob, solver_parameters= params)

sol_it = Function(W, name='sol_it')

def monitor(ksp, iteration_number, norm0):
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

