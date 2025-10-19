from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

parser = ArgumentParser(
    description='Shifted simplified steady Linear Boussinesq equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--nx', type=int, default=40, help='Number of cells along horizontal direction.')
parser.add_argument('--nz', type=int, default=20, help='Number of layers to extrude.')
parser.add_argument('--length', type=float, default=1.0, help='Horizontal length of our solution domain.')
parser.add_argument('--height', type=float, default=1.0, help='Height of our solution domain.')
parser.add_argument('--dt', type=float, default=1.0, help='Time stepping parameter.')
parser.add_argument('--shift', type=float, default=1.0, help='Shift parameter for the shift preconditioner.')
parser.add_argument('--direct', action='store_true', help='If the direct PETSc solver is used to solve the Schur complement.')
parser.add_argument('--shiftpc', action='store_true', help='If True, the equation is preconditioned using the shift preconditioner, if False, the shifted equation is solved.')
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')

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

def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )

class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, u, v):
        print("HDivSchurPC.form() is being called!")
        appctx_PC = self.get_appctx(pc)
        dtc = appctx_PC["dt"]
        delta = appctx_PC["shift"]
        W = u.function_space()
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        w = vector_3D(wxz, wy)
        p = - Constant(1.0) / delta * div(velo)
        Jp = lhs(utils.SLB_velocity(velo, p, b, w, dtc, twoD=False))
        Jp += lhs(utils.SLB_buoyancy(velo, b, q, dtc, twoD=False))
        # k = as_vector([0., 0., 1.])
        # Jp = (inner(velo, w) + dt / shift * div(velo) * div(w))*dx # TODO: The delta shifting parameter enters here.
        # Jp -= dt * inner(w, k) * b * dx
        # Jp += q * b *dx + dt * q * inner(k, velo) * dx
        #  Boundary conditions
        bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
        bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
        bcs = [bc1, bc2]
        return (Jp, bcs)


distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 2)}
deg = 2
m = PeriodicIntervalMesh(nx, length,distribution_parameters=distribution_parameters)
mesh_old = ExtrudedMesh(m, nz, layer_height=height/nz)
x, z = SpatialCoordinate(mesh_old)
coord_fs = VectorFunctionSpace(mesh_old, "DG", deg-1, dim=3)
new_coord = assemble(interpolate(as_vector([x, 0, z]), coord_fs))
mesh = Mesh(new_coord)
mesh.init_cell_orientations(utils.j())
x, y, z = SpatialCoordinate(mesh)

V_2D = utils.extrude_RT(mesh, k=deg) # RT2
DG_km1 = FunctionSpace(mesh, 'DG', deg-1) # DG1
Vy = DG_km1
Pressure = DG_km1
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure

U = Function(W)

uxz, uy, b, p = split(U)
wxz, wy, q, phi = TestFunctions(W)
u = vector_3D(uxz, uy)
w = vector_3D(wxz, wy)

f = Function(V_2D).project(as_vector([sin(2 * pi * x), 0.0, sin(2* pi * z / height)]))
g = Function(DG_km1).interpolate(sin(2*pi*x) + sin(2*pi*z))

# DiricheletBC
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]

eqn = utils.SLB_velocity(u, p, b, w, dt, twoD=False)
eqn -= inner(w, f) * dx
eqn += utils.SLB_buoyancy(u, b, q, dt, twoD=False)
eqn -= g * q * dx
eqn += utils.SLB_pressure(u, phi)
shift_eqn = eqn + shift * p * phi * dx
Jp = derivative(shift_eqn, U)

# Pressure Nullspace
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])

if args.direct:
    params_schur= {
            'mat_type': 'aij',
            'snes_type':'ksponly', 
            'ksp_type': 'gmres', 
            'ksp_view': ':SLB2D.txt',
            'snes_monitor':None, 
            'ksp_monitor':None, 
            'ksp_error_if_not_converged': None,
            'pc_type':'fieldsplit',
            'pc_fieldsplit_type': 'schur',
            'pc_fieldsplit_schur_fact_type': 'full',
            'pc_fieldsplit_schur_precondition': 'full',
            'pc_fieldsplit_0_fields': '3',
            'pc_fieldsplit_1_fields': '0,1,2',
            'fieldsplit_0': {
                'ksp_type': 'preonly',
                'pc_type': 'bjacobi',
                'sub_pc_type': 'ilu',
                # 'ksp_type': 'preonly',
                # 'pc_type': 'lu',
                # 'pc_factor_mat_solver_type': 'petsc',
            },
            'fieldsplit_1': {
                'ksp_type': 'preonly',
                'ksp_monitor': None,
                'mat_view':':mat2.txt',
                # 'pc_type': 'lu',
                # 'pc_factor_mat_solver_type': 'petsc',
                'pc_type':'ksp',
                'ksp_ksp_type': 'preonly',
                'ksp_pc_type':'lu',
                'ksp_ksp_monitor': None,
                'ksp_pc_factor_mat_solver_type':'petsc',
                },
    }
    print("==========================The Schur complement is provided using PETSc.==============")

else:
    pc_params = {
        'mat_view':':mat.txt',
        'pc_type':'ksp',
        'ksp_ksp_type': 'preonly',
        'ksp_pc_type':'lu',
        'ksp_ksp_monitor': None,
        # 'pc_type': 'lu',
        'ksp_pc_factor_mat_solver_type': 'petsc',
    }

    params_schur = {
        'mat_type': 'aij',
        'ksp_type': 'gmres',
        'snes_type':'ksponly',
        # 'ksp_atol': 0,
        # 'ksp_rtol': 1e-8,
        'ksp_view': ':SLB2D.txt',
        'snes_monitor': None,
        'ksp_monitor_true_residual': None,
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        'pc_fieldsplit_schur_precondition': 'full',
        'pc_fieldsplit_0_fields': '3',
        'pc_fieldsplit_1_fields': '0,1,2',
        'fieldsplit_0': { # Doing a pure mass solve for the pressure block.
            'ksp_type': 'preonly',
            'pc_type': 'bjacobi',
            'sub_pc_type': 'ilu',
            # 'pc_factor_mat_solver_type': 'mumps',
        },
        'fieldsplit_1': {
            'ksp_type': 'preonly',
            'ksp_monitor': None,
            'pc_type': 'python',
            'pc_python_type': __name__ + '.HDivSchurPC',
            'helmholtzschurpc': pc_params,
            },
    }
    print("==========================The Schur complement is provided using AuxPC.==============")

if args.shiftpc:
    nprob = NonlinearVariationalProblem(eqn, U, bcs=bcs, Jp=Jp)
    print("================Solving the original equation using the shifted preconditioner.=====================")
else:
    nprob = NonlinearVariationalProblem(shift_eqn, U, bcs=bcs)
    print("================Solving the shifted equation.=====================")
nsolver = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=params_schur, appctx=appctx)

nsolver.solve()
name = 'SLB_slice'
file_slb = VTKFile(name+'.pvd')
uxz, uy, b, p = U.subfunctions
uxz.rename("in-plane-vel")
uy.rename("y-vel")
b.rename("buoyancy")
p.rename("pressure")
f.rename("f_RHS")
g.rename("g_RHS")
file_slb.write(uxz, uy, b, p, f, g)
