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

# ! Parameter settings
parser.add_argument('--nx', type=int, default=40, help='Number of cells along horizontal direction.')
parser.add_argument('--nz', type=int, default=20, help='Number of layers to extrude.')
parser.add_argument('--length', type=float, default=1.0, help='Horizontal length of our solution domain.')
parser.add_argument('--height', type=float, default=1.0, help='Height of our solution domain.')
parser.add_argument('--refinement', type=int, default=2, help='Levels of the multigrid.')
parser.add_argument('--degree', type=int, default=2, help='Order of the element.')
parser.add_argument('--dt', type=float, default=1.0, help='Time stepping parameter.')
parser.add_argument('--tmax', type=float, default=2.0, help='Time period that we solve.')
parser.add_argument('--shift', type=float, default=1.0, help='Shift parameter for the shift preconditioner.')

# ! Test settings
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')
parser.add_argument('--no_rotation', action='store_false', help='If true, no Coriolis term will be imposed in the equation.')
parser.add_argument('--dt_test', action='store_true', help='If true, save the error data storing dt parameters.')
parser.add_argument('--ar_test', action='store_true', help='If true, save the error data storing AR parameters.')
parser.add_argument('--dx_test', action='store_true', help='If true, save the error data storing dx parameters.')
parser.add_argument('--dz_test', action='store_true', help='If true, save the error data storing dz parameters.')
parser.add_argument('--rtol', type=float, default=1.0e-7, help='Relative tolerance for the ksp of linear solver.')
parser.add_argument('--maxit', type=int, default=50, help='Max iteration number for the first ksp of the linear solve.')

# ! Solver settings 
parser.add_argument('--U_mean', type=float, default=0.0, help='Constant horizontal mean flow in x-direction (m/s). Set to 0 for no mean flow.')
parser.add_argument('--direct', action='store_true', help='If true, solve the Schur complement using direct LU.')
parser.add_argument('--timing', action='store_true', help='If true, run the code without monitoring and test for the time.')
parser.add_argument('--reordering', action='store_true', help='If true, run the code with RCM reordering.')
parser.add_argument('--richardson', action='store_true', help='If true, run the code with Richardson iteration for the fieldsplit_1 solve for Schur complement.')
parser.add_argument('--monolithic', action='store_true', help='If true, use a monolithic geometric-multigrid solver with an (extruded) star line patch smoother acting on the whole shifted system, instead of the fieldsplit Schur-complement approach. Takes precedence over --direct.')
parser.add_argument('--cheb_it', type=int, default=6, help='Number of Chebyshev smoother iterations per MG level for the monolithic patch smoother.')

args = parser.parse_known_args()
args = args[0]

if args.show_args:
    PETSc.Sys.Print(args)

use_rotation = not args.no_rotation
monitor_run = not args.timing
any_test = args.dt_test or args.ar_test or args.dx_test or args.dz_test
if args.monolithic:
    solver_name = 'monolithic MG ASMExtrudedStar (star line patch)'
elif args.direct:
    solver_name = 'direct'
else:
    solver_name = 'MG ASMStar'

# if args.timing:
#     opts = PETSc.Options()
#     opts["log_view"] = ":log_view.txt"
#     # opts["log_view_memory"] = None

dt = Constant(args.dt)
tmax = args.tmax
shift = Constant(args.shift)
nx=args.nx
nz=args.nz
length=args.length
height=args.height
deg = args.degree
ar = height / length
deltax = length / nx
deltaz = height / nz
appctx = {
    "dt": dt,
    "shift":shift,
}

print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print(f"Physical domain with length{length} and height {height}, aspect ratio {ar}.")
print(f"Number of elements in x direction {nx} and z direction {nz}.")
print(f"Time stepping parameters dt {args.dt}, total time {tmax}")
print(f"Shift parameter {args.shift}, proportional constant C1 {args.shift * args.dt**1.5}.")
print(f"The code is running with {solver_name} solver for the Schur complement created.")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )

class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, v, u):
        prefix = (pc.getOptionsPrefix() or "") + self._prefix
        rotation = PETSc.Options().getBool(f"{prefix}use_rotation", False) # ! Use this petsc option can make the equation to be consistent.
        appctx_PC = self.get_appctx(pc)
        dtc = appctx_PC["dt"]
        delta = appctx_PC["shift"]
        n = appctx_PC["n"]
        W = u.function_space()
        One = as_vector([1., 1., 1.])
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        unph = Constant(0.5) * (velo + One)
        bnph = Constant(0.5) * (b + Constant(1.))
        w = vector_3D(wxz, wy)
        pnp1 = - Constant(1.) / delta * div(velo)
        Jp = lhs(utils.LB_velocity(velo, One, unph, w, bnph, pnp1, n, dtc, use_rotation=rotation, twoD=False, U_mean=Constant(args.U_mean)))
        Jp += lhs(utils.LB_buoyancy(b, Constant(1.), q, unph, bnph, n, dtc, twoD=False, U_mean=Constant(args.U_mean)))
        #  Boundary conditions
        _, bcs = super().form(pc, u, v)
        return (Jp, bcs)

class MonolithicShiftPC(AuxiliaryOperatorPC):
    """Provide the full shifted Jacobian as an auxiliary operator.

    Wrapping the monolithic operator here (rather than passing it as the
    problem's ``Jp``) gives the inner solver a context with ``Jp = None``
    (J == P), which is what geometric multigrid needs in order to coarsen
    the operator. Applying ``pc_type: mg`` to a global ``Jp`` directly trips
    the ``P.handle == ctx._pjac.petscmat.handle`` assertion in
    firedrake/solving_utils.py during ``PCSetUp_MG``.
    """
    _prefix = "monoshiftpc_"
    def form(self, pc, v, u):
        # Reuse the globally-built shifted Jacobian form, retargeted onto this
        # PC's own test (v) and trial (u) functions.
        test, trial = Jp.arguments()
        a = replace(Jp, {test: v, trial: u})
        _, bcs = super().form(pc, u, v)
        return (a, bcs)

distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length,distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=args.refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (args.refinement+1), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]
finest_mesh_name = "finest"
mesh.name = finest_mesh_name

x, y, z = SpatialCoordinate(mesh)
# ! Lowest order element here. i.e. order 1. with degree = 1. Test with degree = 2 in general.
V_2D = utils.extrude_RT(mesh, k=deg) # ! In common language, RT0 with k=1
Vy = FunctionSpace(mesh, 'DG', deg-1)
Pressure = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure

Un = Function(W)
Unp1 = Function(W)

unxz, uny, bn, pn = split(Un) # ! split for writing the equations.
unp1xz, unp1y, bnp1, pnp1 = split(Unp1)
w_xz, wy, q, phi = TestFunctions(W)

# Initial condition
xc = Constant(length/2)
yc = Constant(length/2)
a = Constant(5000)
U = Constant(0.)
# This is a 4 components function.
u0_slice, u0yic, b0ic, p0ic = Un.subfunctions # ! subfunction for data assignment
u1_slice, u1yic, b1ic, p1ic = Unp1.subfunctions
b0ic.project(3e-4*sin(pi*z/height)/(1+((x-xc)**2)/a**2))
b1ic.project(3e-4*sin(pi*z/height)/(1+((x-xc)**2)/a**2))
# b0ic.project(1e-2*sin(pi*z/height)/(1+((x-xc)**2)/a**2))
# b1ic.project(1e-2*sin(pi*z/height)/(1+((x-xc)**2)/a**2))
# print('===============================================')
# print('Initial condition has been interpolated')
# name = 'ic'
# file_lb = VTKFile(name+'.pvd')
# u0, u0y, b0, P0 = Un.subfunctions
# file_lb.write(u0, u0y, b0, P0)
# print("Save initial condition.")

# DiricheletBC
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]

un = vector_3D(unxz, uny)
unp1 = vector_3D(unp1xz, unp1y)
unph = 0.5 * (un+unp1)
bnph = 0.5 * (bn+bnp1)
pnph = 0.5 * (pn+pnp1)
w = vector_3D(w_xz, wy)
n = FacetNormal(mesh)
appctx.update({"n": n})

eqn = utils.LB_velocity(unp1, un, unph, w, bnph, pnp1, n, dt, use_rotation=use_rotation, U_mean=Constant(args.U_mean))
eqn += utils.LB_buoyancy(bnp1, bn, q, unph, bnph, n, dt, U_mean=Constant(args.U_mean))
eqn += utils.LB_pressure(unp1, phi)
shift_eqn = eqn + shift * pnp1 * phi * dx
Jp = derivative(shift_eqn, Unp1)

# Pressure Nullspace
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])

if args.direct:
    helmholtz_schur_pc_params = {
        'pc_type':'ksp',
        'ksp_ksp_type': 'preonly',
        # 'mat_view':':matSchurAuxPC.txt',
        'ksp_pc_type':'lu',
        # 'ksp_ksp_monitor': None,
    }

else:
    helmholtz_schur_pc_params = {
            # 'ksp_type': 'preonly',
            # 'ksp_max_its': 30,
            'pc_type': 'mg',
            'pc_mg_type': 'full',
            'pc_mg_cycle_type':'v',
            'mg_levels': {
                # 'ksp_type': 'gmres',
                # 'ksp_type':'richardson',
                # 'ksp_type': 'chebyshev',
                # 'ksp_richardson_scale': 0.5,
                # 'ksp_richardson_self_scale':None,
                # 'ksp_max_it': 1, # ? more robust for larger max_it here.
                # 'ksp_monitor':None,
                "pc_type": "python",
                "pc_python_type": "firedrake.ASMStarPC",
                "pc_star_construct_dim": 0,
                "pc_star_sub_sub_pc_type": "lu",
                "pc_star_view_patch_sizes":True,
                # 'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
                # 'pc_star_sub_sub_pc_factor_reuse_ordering': None,
                # 'pc_star_sub_sub_pc_factor_mat_solver_type': 'mumps',
                # 'pc_star_sub_sub_pc_factor_mat_solver_type': 'superlu_dist',
                # "pc_star_sub_sub_pc_type": "svd",
                # "pc_star_sub_sub_pc_svd_monitor": None,
            },
            'mg_coarse': {
                'ksp_type': 'preonly',
                'pc_type': 'lu',
            },
        }
    if args.reordering:
        helmholtz_schur_pc_params.update({
            'mg_levels_pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
            'mg_levels_pc_star_sub_sub_pc_factor_reuse_ordering': None,
        })
    if args.richardson:
        helmholtz_schur_pc_params.update({
            # 'mg_levels_ksp_type': 'richardson',
            'mg_levels_ksp_type':'chebyshev', # ! chebyshev smoother is working.
            # 'mg_levels_ksp_type': 'gmres',
            'mg_levels_ksp_max_it':6,
        })
    else:
        helmholtz_schur_pc_params.update({
            'mg_levels_ksp_type': 'gmres',
            'mg_levels_ksp_max_it':6,
        })

params_schur = {
    # 'mat_type': 'aij',
    # 'log_view':':log_view.txt',
    # 'log_view_memory':':log_view_memory.txt',

    # 'ksp_type': 'fgmres', # ! this can also be tuned.
    'snes_type':'ksponly',
    'ksp_atol': 0,
    'ksp_rtol': args.rtol,
    'ksp_max_it': args.maxit,
    'ksp_converged_maxits': None, # ! When max_it is reached, setting this will pass the convergence test and make the solver run, instead of raising a ConvergenceError. Distinguish the type of convergence in ConvergedReason instead!
    'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_fact_type': 'full',
    'pc_fieldsplit_0_fields': '3',
    'pc_fieldsplit_1_fields': '0,1,2',
    'fieldsplit_0': { # Doing a pure mass solve for the pressure block.
        'ksp_type': 'preonly',
        'pc_type': 'bjacobi',
        'sub_pc_type': 'ilu',
        # 'pc_factor_mat_solver_type': 'mumps',
    },
    'fieldsplit_1': {
        'use_rotation':use_rotation,
        # 'ksp_type': 'fgmres', # ! need to tune this.
        # 'ksp_type': 'richardson',
        # 'ksp_richardson_scale': 1.0,
        # 'ksp_richardson_self_scale':None,
        # 'ksp_atol': 0,
        # 'ksp_rtol': 1e-7, # ? Do I need to set this?
        # 'mat_view':':field_1_mat_aux.txt',
        'pc_type': 'python',
        'pc_python_type': __name__ + '.HDivSchurPC',
        'helmholtzschurpc': helmholtz_schur_pc_params,
    },
}
params_schur['ksp_type'] = 'gmres' if args.richardson else 'fgmres'
params_schur['fieldsplit_1_ksp_type'] = 'preonly' if args.richardson else 'fgmres'

# ! Monolithic solver: geometric multigrid applied to the *whole* shifted system Jp
# ! (supplied via MonolithicShiftPC), with a star patch smoother, optional RCM
# ! reordering on the patch sub-solves and a Chebyshev smoother KSP wrapping the patch.
# ! Use firedrake.ASMExtrudedStarPC below for vertical line/column patches instead.

mg_levels_mono = {
    # 'ksp_type': 'chebyshev',
    'ksp_type':'gmres',     # Chebyshev iteration on the patch smoother.
    'ksp_max_it': 6,
    'pc_type': 'python',
    'pc_python_type': 'firedrake.ASMStarPC',  # switch to firedrake.ASMExtrudedStarPC for vertical line patches
    'pc_star_construct_dim': 0,
    'pc_star_sub_sub_pc_type': 'lu',
    'pc_star_view_patch_sizes': True,
}
if args.reordering:
    mg_levels_mono.update({
        'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
        'pc_star_sub_sub_pc_factor_reuse_ordering': None,
    })

params_monolithic = {
    'snes_type': 'ksponly',
    'ksp_type': 'fgmres',
    'ksp_atol': 0,
    'ksp_rtol': args.rtol,
    'ksp_max_it': args.maxit,
    'ksp_converged_maxits': None,
    # Shifted operator supplied via MonolithicShiftPC so the inner MG context
    # has J == P, which geometric multigrid requires to coarsen the operator.
    'pc_type': 'python',
    'pc_python_type': __name__ + '.MonolithicShiftPC',
    'monoshiftpc': {
        'mat_type': 'aij',  # patch PCs need an assembled (non-nested) operator.
        'pc_type': 'mg',
        'pc_mg_type': 'full',
        'pc_mg_cycle_type': 'v',
        'mg_levels': mg_levels_mono,
        'mg_coarse': {
            'ksp_type': 'preonly',
            'pc_type': 'lu',
        },
    },
}
# ! Choose which solver configuration to run with. Monolithic takes precedence.
params = params_monolithic if args.monolithic else params_schur

if args.timing:
    params['ksp_view'] = ':slice3D.txt'
else:
    params.update({
        'snes_monitor': None,
        # 'ksp_monitor': None,
        'ksp_converged_rate': None,
        'ksp_monitor_true_residual': None,
        # "ksp_error_if_not_converged": False,
        # "snes_error_if_not_converged": False,
    })
    if not args.monolithic:
        params['fieldsplit_1'].update({
            'ksp_monitor': None,
            'ksp_converged_reason': f':fieldsplit1_ksp_dt{args.dt}_shift{args.shift}.txt',
        })

# 'fieldsplit_1_ksp_type': 'richardson', 'fieldsplit_1_ksp_richardson_scale': 1.0,
# print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
nprob = NonlinearVariationalProblem(eqn, Unp1, bcs=bcs, Jp=Jp)
# nprob = NonlinearVariationalProblem(shift_eqn, Unp1, bcs=bcs) # this will set the non-shifted equation.
nsolver = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=params, appctx=appctx)

if monitor_run:
    # Set checkpointing for saving the data.
    error_list = []
    sol_it = Function(W, name='sol_it')
    sol_final = Function(W, name='sol_final')
    # with CheckpointFile('sol_mesh.h5', 'w') as chk:
    #     chk.save_mesh(mesh)

    def monitor(ksp, iteration_number, norm0):
        # print('The monitor starts to build the solution.')
        sol = ksp.buildSolution()
        with sol_it.dat.vec_wo as it_vec:
            sol.copy(result=it_vec)
        # print('The monitor starts to calculate the error.')
        error = norm(sol_it - sol_final) / norm(sol_final)
        error_list.append(error)

# Time Stepping
if monitor_run:
    name = 'lb_slice_imp_ASM'
    file_lb = VTKFile(name+'.pvd')
    un, uny, bn, pn = Un.subfunctions
    un.rename("in-plane-vel")
    uny.rename("y-vel")
    bn.rename("buoyancy")
    pn.rename("pressure")
    file_lb.write(un, uny, bn, pn)
Unp1.assign(Un)
t = 0.0
dumpt = args.dt
tdump = 0.
j = 0
while t < tmax - 0.5 * args.dt:
    print(f"=======================================The solver is currently solving for time:{t}==========================")
    t += args.dt
    tdump += args.dt
    if j == 0:
        nsolver.solve()
    else:
        U_restart = Unp1.copy(deepcopy=True)
        with PETSc.Log.Stage("Warm-up-solve"):
            nsolver.solve()
        if not monitor_run: # ! test the CPU Time.
            Unp1.assign(U_restart) # ! Assign the original velocity to restart the solver.
            with PETSc.Log.Stage("Official-Run"):
                nsolver.solve()
        if monitor_run and any_test:
            final_sol = nsolver.snes.ksp.buildSolution()
            with sol_final.dat.vec_wo as final_vec:
                final_sol.copy(result=final_vec)
            Unp1.assign(U_restart) # ! Assign the original velocity to restart the solver.
            nsolver.snes.ksp.setMonitor(monitor)
            with PETSc.Log.Stage("Official-Run"):
                nsolver.solve()
            reason = nsolver.snes.ksp.getConvergedReason()
            print("*************************************************", reason)
            converged_it_num = nsolver.snes.ksp.getIterationNumber()
            print(f"KSP is converged in {converged_it_num} iterations and monitor is working on time step {j}.")
            if args.dt_test:
                np.savetxt(f'error_dt{args.dt}_shift{args.shift}.out', error_list)
            if args.ar_test:
                np.savetxt(f'error_dt{args.dt}_ar{ar}.out', error_list)
            if args.dx_test:
                np.savetxt(f'error_dt{args.dt}_nx{nx}.out', error_list)
            if args.dz_test:
                np.savetxt(f'error_dt{args.dt}_nz{nz}.out', error_list)
    Un.assign(Unp1)
    j += 1
    if tdump > dumpt - args.dt*0.5:
        if monitor_run:
            file_lb.write(un, uny, bn, pn)
        tdump -= dumpt
