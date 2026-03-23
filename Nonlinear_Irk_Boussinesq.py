from firedrake import *
from irksome import GaussLegendre, Dt, MeshConstant, TimeStepper, IRKAuxiliaryOperatorPC
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter
from firedrake.dmhooks import get_appctx as get_snesctx

parser = ArgumentParser(
    description='Shifted simplified steady Linear Boussinesq equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--nx', type=int, default=40, help='Number of cells along horizontal direction.')
parser.add_argument('--nz', type=int, default=20, help='Number of layers to extrude.')
parser.add_argument('--length', type=float, default=1.0, help='Horizontal length of our solution domain.')
parser.add_argument('--height', type=float, default=1.0, help='Height of our solution domain.')
parser.add_argument('--refinement', type=int, default=2, help='Levels of the multigrid.')
parser.add_argument('--degree', type=int, default=1, help='Order of the element.')
parser.add_argument('--dt', type=float, default=1.0, help='Time stepping parameter.')
parser.add_argument('--tmax', type=float, default=2.0, help='Time period that we solve.')
parser.add_argument('--shift', type=float, default=1.0, help='Shift parameter for the shift preconditioner.')
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')
parser.add_argument('--no_rotation', action='store_false', help='If true, no Coriolis term will be imposed in the equation.')
parser.add_argument('--rtol', type=float, default=1.0e-9, help='Relative tolerance for the ksp of linear solver.')
parser.add_argument('--maxit', type=int, default=150, help='Max iteration number for the first ksp of the linear solve.')

args = parser.parse_known_args()
args = args[0]

if args.show_args:
    PETSc.Sys.Print(args)

if args.no_rotation:
    use_rotation = False
else:
    use_rotation = True


solver_name = 'MG ASMStar'

nx=args.nx
nz=args.nz
length=args.length
height=args.height
deg = args.degree
ar = height / length
deltax = length / nx
deltaz = height / nz
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print(f"Physical domain with length{length} and height {height}, aspect ratio {ar}.")
print(f"Number of elements in x direction {nx} and z direction {nz}.")
print(f"Time stepping parameters dt {args.dt}, total time {args.tmax}")
print(f"Shift parameter {args.shift}, proportional constant C1 {args.shift * args.dt**1.5}.")
print(f"The code is running with {solver_name} solver for the Schur complement created.")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )

distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length,distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=args.refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (args.refinement+1), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3) 
mesh = new_mh[-1]
finest_mesh_name = "finest"
mesh.name = finest_mesh_name

MC = MeshConstant(mesh)
dt = MC.Constant(args.dt)
shift = Constant(args.shift)
tmax = args.tmax
t = MC.Constant(0.0)
appctx = {
    "dt": dt,
    "shift":shift,
}

x, y, z = SpatialCoordinate(mesh)
n = FacetNormal(mesh)
appctx.update({"n": n})
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg-1)
Pressure = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure

U = Function(W)

uxz, uy, b, p = split(U) # ! split for writing the equations.
w_xz, wy, q, phi = TestFunctions(W)

# Initial condition
xc = Constant(length/2)
yc = Constant(length/2)
a = Constant(5000)
# U_mean = Constant(0.)
# This is a 4 components function.
u0_slice, u0yic, b0ic, p0ic = U.subfunctions # ! subfunction for data assignment
b0ic.project(sin(pi*z/height)/(1+((x-xc)**2)/a**2))

# DiricheletBC
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]

u = vector_3D(uxz, uy)
w = vector_3D(w_xz, wy)

eqn = utils.Nonlinear_velocity_Irk(u, w, b, p, dt, n, use_rotation=use_rotation)
eqn += utils.Nonlinear_buoyancy_Irk(b, q, u, dt, n)
eqn += utils.Nonlinear_pressure_Irk(u, phi)

# Pressure Nullspace
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])

class HDivSchurPC(IRKAuxiliaryOperatorPC):
    _prefix = 'shiftedschurpc_'
    def getNewForm(self, pc, u0, test):
        prefix = (pc.getOptionsPrefix() or "") + self._prefix
        rotation = PETSc.Options().getBool(f"{prefix}use_rotation", False) # ! Use this petsc option can make the equation to be consistent.
        appctx_PC = self.get_appctx(pc) # ! returning dmhooks.get_appctx(pc).appctx.
        dtc = appctx_PC["dt"]
        delta = appctx_PC["shift"]
        n = appctx_PC["n"]
        W = u0.function_space()
        uxz, uy, b, p = split(u0)
        wxz, wy, q, phi = split(test)
        u = vector_3D(uxz, uy)
        w = vector_3D(wxz, wy)
        # ? The pressure elimination happened here, and no more pressure equation.
        # p = - Constant(1.) / delta * div(u)
        
        F = utils.Nonlinear_velocity_Irk(u, w, b, p, dtc, n, use_rotation=rotation)
        F += utils.Nonlinear_buoyancy_Irk(b, q, u, dtc, n)
        F += utils.Nonlinear_pressure_Irk(u, phi)
        F += delta * p * phi * dx
        # print("::::::::::::::::::::::")
        #  Boundary conditions
        #  Boundary conditions
        bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
        bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
        bcs = [bc1, bc2]
        # _, bcs = super().form(pc, u0, test)
        return (F, bcs)




# helmholtz_schur_pc_params = {
#         # 'ksp_type': 'preonly',
#         # 'ksp_max_its': 30,
#         'pc_type': 'mg',
#         'pc_mg_type': 'full',
#         'pc_mg_cycle_type':'v',
#         'mg_levels': {
#             'ksp_type': 'gmres',
#             'ksp_max_it': 6, # ? more robust for larger max_it here.
#             # 'ksp_monitor':None,
#             "pc_type": "python",
#             "pc_python_type": "firedrake.ASMStarPC",
#             "pc_star_construct_dim": 0,
#             "pc_star_sub_sub_pc_type": "lu",
#             'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
#             'pc_star_sub_sub_pc_factor_reuse_ordering': None,
#         },
#         'mg_coarse': {
#             'ksp_type': 'preonly',
#             'pc_type': 'lu',
#         },
#     }

# ? Using the fieldsplit Schur complement to solve the preconditioned Schur complement and use the direct solver to solve the Schur complement.
# shifted_schur_pc_params ={
#     'pc_type':'ksp',
#     'ksp_ksp_type': 'preonly',
#     'ksp_pc_type':'lu',
#     'ksp_ksp_monitor': None,
#     # "pc_type": "lu",
#     # "pc_factor_mat_solver_type": "mumps",
# }


# params_schur = {
#     'mat_type': 'matfree',
#     'ksp_view': ':Nonlinear_slice3D.txt',
#     'ksp_type': 'gmres',
#     'snes_type':'ksponly',
#     'ksp_atol': 0,
#     'ksp_rtol': args.rtol,
#     'ksp_max_it': args.maxit,
#     'ksp_converged_maxits': None,
#     'snes_monitor': None,
#     # 'ksp_monitor': None,
#     'ksp_converged_rate':None,
#     'ksp_monitor_true_residual': None,
#     # "ksp_error_if_not_converged": False,
#     # "snes_error_if_not_converged": False,
#     'pc_type': 'fieldsplit',
#     'pc_fieldsplit_type': 'schur',
#     'pc_fieldsplit_schur_fact_type': 'full',
#     'pc_fielsplit_schur_precondition': 'a11',
#     'pc_fieldsplit_0_fields': '3',
#     'pc_fieldsplit_1_fields': '0,1,2',
#     'fieldsplit_0': { # Doing a pure mass solve for the pressure block.
#         'ksp_type': 'preonly',
#         'pc_type':'python',
#         'pc_python_type':'firedrake.AssembledPC',
#         'assembled_pc_type': 'bjacobi',
#         'assembled_sub_pc_type': 'ilu', # ILU needs an assembled matrix so that AssembledPC is needed.
#     },
#     'fieldsplit_1': {
#         'helmholtzschurpc_use_rotation':use_rotation,
#         'ksp_type': 'preonly', # ! need to tune this.
#         'ksp_monitor': None,
#         'ksp_converged_reason': f':fieldsplit1_ksp_dt{args.dt}_shift{args.shift}.txt',
#         # 'ksp_atol': 0,
#         # 'ksp_rtol': 1e-7, # ? Do I need to set this?
#         # 'mat_view':':field_1_mat_aux.txt',
#         'pc_type': 'python',
#         'pc_python_type': __name__ + '.HDivSchurPC',
#         'shiftedschurpc': shifted_schur_pc_params,
#         },
# }


# Direct solve for only the shifted equation using the AuxOPPC.

shifted_schur_pc_params ={
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
}
params_schur = {
    'mat_type': 'matfree', # ! Auxiliary Operator PC needs a mat-free operator.
    'ksp_view': ':Nonlinear_slice3D.txt',

    'ksp_type': 'gmres',
    'snes_type':'newtonls', # Newton Line Search
    'ksp_atol': 0,
    'ksp_rtol': args.rtol,
    'ksp_max_it': args.maxit,
    'ksp_converged_maxits': None,
    'snes_monitor': None,
    # 'ksp_monitor': None,
    'ksp_converged_rate':None,
    'ksp_monitor_true_residual': None,
    # "ksp_error_if_not_converged": False,
    # "snes_error_if_not_converged": False,
    'pc_type': 'python',
    'pc_python_type': __name__ + '.HDivSchurPC',
    'shiftedschurpc': shifted_schur_pc_params,
}


butcher_tableau = GaussLegendre(1)

stepper = TimeStepper(eqn, butcher_tableau, t, dt, U, bcs=bcs, solver_parameters=params_schur, appctx=appctx)


name = 'Nonlinear_slice_imp_ASM'
file_lb = VTKFile(name+'.pvd')
un, uny, bn, pn = U.subfunctions
un.rename("in-plane-vel")
uny.rename("y-vel")
bn.rename("buoyancy")
pn.rename("pressure")
file_lb.write(un, uny, bn, pn)

dumpt = args.dt
tdump = 0.
while (float(t) < tmax - 0.5 * args.dt):
    tdump += args.dt
    stepper.advance()
    print(float(t))
    t.assign(float(t) + float(dt))
    if tdump > dumpt - 0.5*args.dt:
        file_lb.write(un, uny, bn, pn)
    tdump -= dumpt



