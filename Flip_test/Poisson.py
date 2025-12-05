from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

parser = ArgumentParser(
    description= 'Poisson equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--nx', type=int, default=40, help='Number of cells along horizontal direction.')
parser.add_argument('--nz', type=int, default=20, help='Number of layers to extrude.')
parser.add_argument('--length', type=float, default=1.0, help='Horizontal length of our solution domain.')
parser.add_argument('--height', type=float, default=1.0, help='Height of our solution domain.')
parser.add_argument('--dt', type=float, default=1.0, help='Time stepping parameter.')
parser.add_argument('--shift', type=float, default=1.0, help='Shift parameter for the shift preconditioner.')
parser.add_argument('--direct', action='store_true', help='If True, the direct PETSc solver is used to solve the Schur complement.')
parser.add_argument('--shiftpc', action='store_true', help='If True, the Poisson equation is preconditioned using the shift preconditioner, if False, the shifted equation is solved.')
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')
parser.add_argument('--threeD', action='store_true', help='If True, the 3D slice problem is solved using embedded 3D mesh.')

args = parser.parse_known_args()
args = args[0]

if args.show_args:
    PETSc.Sys.Print(args)

dt = Constant(args.dt)
shift = Constant(args.shift)
nx=args.nx
nz=args.nz
length=args.length
height=args.height
appctx = {
    "dt": dt,
    "shift":shift,
}
if args.threeD:
    dimension = 3
    print(f"=================The problem is in {dimension}D.")
else:
    dimension = 2
    print(f"=================The problem is in {dimension}D.")

def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )

if args.threeD:
    class HDivSchurPC(AuxiliaryOperatorPC):
        _prefix = "helmholtzschurpc_"
        def form(self, pc, u, v):
            appctx_PC = self.get_appctx(pc)
            dtc = appctx_PC["dt"]
            delta = appctx_PC["shift"]
            W = u.function_space()
            uxz, uy = split(u)
            wxz, wy = split(v)
            velo = uxz + uy * utils.j()
            w = wxz + wy * utils.j()
            Jp = inner(velo, w) * dx + dtc / delta * div(velo) * div(w) * dx
            bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
            bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
            bcs = [bc1, bc2]
            return (Jp, bcs)
else:
    class HDivSchurPC(AuxiliaryOperatorPC):
        _prefix = "helmholtzschurpc_"
        def form(self, pc, u, v):
            appctx_PC = self.get_appctx(pc)
            dtc = appctx["dt"]
            delta = appctx["shift"]
            W = u.function_space()
            Jp = inner(u, v) * dx + dtc / delta * div(u) * div(v) * dx
            bc1 = DirichletBC(W.sub(0), as_vector([0., 0.]), "top")
            bc2 = DirichletBC(W.sub(0), as_vector([0., 0.]), "bottom")
            bcs = [bc1, bc2]
            return (Jp, bcs)


distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 2)}
m = PeriodicIntervalMesh(nx, length,distribution_parameters=distribution_parameters)
deg = 2
if args.threeD:
    mesh_old = ExtrudedMesh(m, nz, layer_height=height/nz)
    x, z = SpatialCoordinate(mesh_old)
    coord_fs = VectorFunctionSpace(mesh_old, "DG", deg-1, dim=3)
    new_coord = assemble(interpolate(as_vector([x, 0, z]), coord_fs))
    mesh = Mesh(new_coord)
    mesh.init_cell_orientations(utils.j())
    x, y, z = SpatialCoordinate(mesh)
else:
    mesh = ExtrudedMesh(m, nz, layer_height=height/nz)
    x, z = SpatialCoordinate(mesh)

V_2D = utils.extrude_RT(mesh, k=deg) # RT2
DG_km1 = FunctionSpace(mesh, 'DG', deg-1) # DG1
Vy = DG_km1
Pressure = DG_km1
if args.threeD:
    W = V_2D * Vy * Pressure
    U = Function(W)
    uxz, uy, p = split(U)
    wxz, wy, phi = TestFunctions(W)
    u = vector_3D(uxz, uy)
    w = vector_3D(wxz, wy)
    F = Function(V_2D).project(as_vector([sin(2 * pi * x), 0.0, sin(2* pi * z / height)]))
    F.rename("f_RHS")
    ic = VTKFile("ic.pvd")
    ic.write(F)
else:
    W = V_2D * Pressure
    U = Function(W)
    u, p = split(U)
    w, phi = TestFunctions(W)
    F = Function(V_2D).project(as_vector([sin(2 * pi * x), sin(2* pi * z / height)]))
    F.rename("f_RHS")
    ic = VTKFile("ic.pvd")
    ic.write(F)

if args.threeD:
    bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
    bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
else:
    bc1 = DirichletBC(W.sub(0), as_vector([0., 0.]), "top")
    bc2 = DirichletBC(W.sub(0), as_vector([0., 0.]), "bottom")
bcs = [bc1, bc2]

eqn = inner(u,w)*dx - dt*div(w)*p*dx - inner(F,w)*dx
eqn += phi*div(u)*dx
shift_eqn = eqn + shift*p*phi*dx
Jp = derivative(shift_eqn, U)

v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
if args.threeD:
    nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), v_basis])
else:
    nullspace = MixedVectorSpaceBasis(W, [W.sub(0), v_basis])


if args.direct:
    params_schur= {
            'mat_type': 'aij',
            'snes_type':'ksponly', 
            'ksp_type': 'gmres', 
            'ksp_view': ':Poisson.txt',
            'snes_monitor':None, 
            'ksp_monitor':None, 
            'ksp_error_if_not_converged': None,
            'pc_type':'fieldsplit',
            'pc_fieldsplit_type': 'schur',
            'pc_fieldsplit_schur_fact_type': 'full',
            'pc_fieldsplit_schur_precondition': 'full',
            'fieldsplit_0': {
                'ksp_type': 'preonly',
                'pc_type': 'bjacobi',
                'sub_pc_type': 'ilu',
            },
            'fieldsplit_1': {
                'mat_view':':matPoissonpetsc.txt',
                'ksp_type': 'preonly',
                'ksp_monitor': None,
                'pc_type': 'lu',
                'pc_factor_mat_solver_type': 'petsc',
                },
    }
    print("==========================The Schur complement is provided using PETSc.==============")

else:
    helmholtz_schur_pc_params = {
        'pc_type':'ksp',
        'ksp_ksp_type': 'preonly',
        'ksp_pc_type':'lu',
        'ksp_ksp_monitor': None,
        # 'pc_type': 'lu',
        # 'pc_factor_mat_solver_type': 'petsc',
    }

    params_schur = {
        'mat_type': 'aij',
        'ksp_type': 'gmres',
        'snes_type':'ksponly',
        # 'ksp_atol': 0,
        # 'ksp_rtol': 1e-8,
        'ksp_view': ':Poisson.txt',
        'snes_monitor': None,
        'ksp_monitor_true_residual': None,
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        'pc_fieldsplit_schur_precondition': 'full',
        'fieldsplit_0': { # Doing a pure mass solve for the pressure block.
            'ksp_type': 'preonly',
            'pc_type': 'bjacobi',
            'sub_pc_type': 'ilu',
            # 'pc_factor_mat_solver_type': 'mumps',
        },
        'fieldsplit_1': {
            'ksp_type': 'preonly',
            'mat_view':':matPoisson.txt',
            'ksp_monitor': None,
            'pc_type': 'python',
            'pc_python_type': __name__ + '.HDivSchurPC',
            'helmholtzschurpc': helmholtz_schur_pc_params,
            },
    }
    print("==========================The Schur complement is provided using AuxPC.==============")

if args.threeD:
    params_schur.update({'pc_fieldsplit_0_fields': '2', 'pc_fieldsplit_1_fields': '0,1',})
else:
    params_schur.update({'pc_fieldsplit_0_fields': '1', 'pc_fieldsplit_1_fields': '0',})


if args.shiftpc:
    nprob = NonlinearVariationalProblem(eqn, U, bcs=bcs, Jp=Jp)
    print("================Solving the original equation using the shifted preconditioner.=====================")
else:
    nprob = NonlinearVariationalProblem(shift_eqn, U, bcs=bcs)
    print("================Solving the shifted equation.=====================")
nsolver = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=params_schur, appctx=appctx)

nsolver.solve()
name = f'Poisson_{dimension}D'
file = VTKFile(name+'.pvd')
if args.threeD:
    uxz, uy, p = U.subfunctions
    uxz.rename("in-plane-vel")
    uy.rename("y-vel")
    p.rename("pressure")
    file.write(uxz, uy, p, F)
else:
    uxz, p = U.subfunctions
    uxz.rename("in-plane-vel")
    p.rename("pressure")
    file.write(uxz, p, F)