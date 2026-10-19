'''
Stand-alone copy of the solver set-up of Nonlinear_Boussinesq_Irk_SC.py (mesh, spaces,
forms from utils.py, HDivSchurPC, solver parameters), used ONLY to validate the reduced
model emulator.py.  It does not modify anything in the repository.

Differences to the production script (all optional, selected on the command line):
  --linear        use the linear forms utils.LB_*_Irk (no advection) instead of the nonlinear
                  ones, i.e. exactly the operator that emulator.py models (with f = 0);
  --inner_lu      replace the multigrid on field 1 by a direct solve (upper bound on what any
                  patch can achieve);
  --smoother_its  number of GMRES iterations of the smoother (6 in the production script).
The second time step is solved twice as in the production script (warm-up, then monitored).

Run from a scratch directory, e.g.
  mpiexec --bind-to none -n 2 python solver_check.py --nx 20 --nz 20 --length 3e5 --height 4e3 \\
      --dt 50 --shift 1.096622711232e-08 --patch vanka --linear
'''
import sys
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
from firedrake import *                                                   # noqa
from irksome import GaussLegendre, Dt, MeshConstant, TimeStepper, IRKAuxiliaryOperatorPC
import numpy as np
from petsc4py import PETSc
import utils
from argparse import ArgumentParser

print = PETSc.Sys.Print

parser = ArgumentParser()
parser.add_argument('--nx', type=int, default=20)
parser.add_argument('--nz', type=int, default=20)
parser.add_argument('--length', type=float, default=3.0e5)
parser.add_argument('--height', type=float, default=4.0e3)
parser.add_argument('--refinement', type=int, default=2)
parser.add_argument('--degree', type=int, default=2)
parser.add_argument('--dt', type=float, default=50.0)
parser.add_argument('--shift', type=float, default=1.096622711232e-08)
parser.add_argument('--rtol', type=float, default=1.0e-6)
parser.add_argument('--atol', type=float, default=1.0e-9)
parser.add_argument('--maxit', type=int, default=100)
parser.add_argument('--patch', type=str, default='star', choices=['star', 'vanka'])
parser.add_argument('--linear', action='store_true')
parser.add_argument('--inner_lu', action='store_true')
parser.add_argument('--smoother_its', type=int, default=6)
parser.add_argument('--amplitude', type=float, default=3.0e-4)
parser.add_argument('--tag', type=str, default='')
args, _ = parser.parse_known_args()

nx, nz, length, height, deg = args.nx, args.nz, args.length, args.height, args.degree


def vector_3D(u, uy):
    return u + uy * utils.j()


distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=args.refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (args.refinement + 1), extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]

MC = MeshConstant(mesh)
dt = MC.Constant(args.dt)
shift = Constant(args.shift)
t = MC.Constant(0.0)
n = FacetNormal(mesh)
appctx = {"dt": dt, "shift": shift, "n": n}

x, y, z = SpatialCoordinate(mesh)
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg - 1)
Pressure = FunctionSpace(mesh, 'DG', deg - 1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure

U = Function(W)
uxz, uy, b, p = split(U)
w_xz, wy, q, phi = TestFunctions(W)
xc = Constant(length / 2)
a = Constant(5000)
u0_slice, u0yic, b0ic, p0ic = U.subfunctions
b0ic.project(args.amplitude * sin(pi * z / height) / (1 + ((x - xc) ** 2) / a ** 2))

bcs = [DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top"),
       DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")]


def forms(u, w, b, p, q, phi, n):
    if args.linear:
        # utils.LB_*_Irk with the Coriolis term removed (f = 0)
        F = inner(w, Dt(u)) * dx - div(w) * p * dx - inner(w, utils.k()) * b * dx
        F += q * Dt(b) * dx + utils.buo_freq() * q * inner(utils.k(), u) * dx
        F += phi * div(u) * dx
        return F
    F = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=False)
    F += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
    F += utils.Nonlinear_pressure_Irk(u, phi)
    return F


u = vector_3D(uxz, uy)
w = vector_3D(w_xz, wy)
eqn = forms(u, w, b, p, q, phi, n)


class HDivSchurPC(IRKAuxiliaryOperatorPC):
    _prefix = 'shiftedschurpc_'

    def getNewForm(self, pc, u0, test):
        appctx_PC = self.get_appctx(pc)
        delta = appctx_PC["shift"]
        n = appctx_PC["n"]
        W = u0.function_space()
        uxz, uy, b, p = split(u0)
        wxz, wy, q, phi = split(test)
        u = vector_3D(uxz, uy)
        w = vector_3D(wxz, wy)
        p = p - Constant(1.) / delta * div(u)
        F = forms(u, w, b, p, q, phi, n)
        F += delta * p * phi * dx
        bcs = [DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top"),
               DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")]
        return (F, bcs)


if args.patch == 'star':
    patch_params = {
        "pc_type": "python",
        "pc_python_type": "firedrake.ASMStarPC",
        "pc_star_construct_dim": 0,
        "pc_star_sub_sub_pc_type": "lu",
        'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
        'pc_star_sub_sub_pc_factor_reuse_ordering': None,
    }
else:
    patch_params = {
        "pc_type": "python",
        "pc_python_type": "firedrake.ASMVankaPC",
        "pc_vanka_construct_dim": 1,
        "pc_vanka_sub_sub_pc_type": "lu",
        'pc_vanka_sub_sub_pc_factor_mat_ordering_type': 'rcm',
        'pc_vanka_sub_sub_pc_factor_reuse_ordering': None,
    }

if args.inner_lu:
    fs1 = {'ksp_type': 'preonly', 'pc_type': 'lu', 'pc_factor_mat_solver_type': 'mumps'}
else:
    fs1 = {
        'ksp_type': 'preonly',
        'pc_type': 'mg',
        'pc_mg_type': 'full',
        'pc_mg_cycle_type': 'v',
        'mg_levels': {'ksp_type': 'gmres', 'ksp_max_it': args.smoother_its, **patch_params},
        'mg_coarse': {'ksp_type': 'preonly', 'pc_type': 'lu'},
    }

shifted_schur_pc_params = {
    'use_rotation': False,
    'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_fact_type': 'full',
    'pc_fieldsplit_schur_precondition': 'a11',
    'pc_fieldsplit_0_fields': '3',
    'pc_fieldsplit_1_fields': '0,1,2',
    'fieldsplit_0': {
        'ksp_type': 'preonly',
        'ksp_reuse_preconditioner': None,
        'pc_type': 'python',
        'pc_python_type': 'firedrake.AssembledPC',
        'assembled_pc_type': 'bjacobi',
        'assembled_sub_pc_type': 'ilu',
    },
    'fieldsplit_1': fs1,
}

params_schur = {
    'mat_type': 'matfree',
    'snes_type': 'newtonls',
    'ksp_type': 'fgmres',
    'ksp_pc_side': 'right',
    'ksp_atol': args.atol,
    'ksp_rtol': args.rtol,
    'snes_max_it': 10,
    'snes_atol': 1e-9,
    'ksp_max_it': args.maxit,
    'ksp_converged_maxits': None,
    'pc_type': 'python',
    'pc_python_type': __name__ + '.HDivSchurPC',
    'shiftedschurpc': shifted_schur_pc_params,
    'snes_monitor': None,
    'ksp_monitor_true_residual': None,
    'ksp_converged_reason': None,
}

stepper = TimeStepper(eqn, GaussLegendre(1), t, dt, U, bcs=bcs, solver_parameters=params_schur, appctx=appctx)

cum = []


def snes_monitor(snes, it, fnorm):
    cum.append((it, fnorm, snes.getLinearSolveIterations()))


print(f'### case: patch={args.patch} linear={args.linear} inner_lu={args.inner_lu} nx={nx} nz={nz} '
      f'L={length} H={height} dt={args.dt} shift={args.shift} smoother_its={args.smoother_its} '
      f'ranks={COMM_WORLD.size} {args.tag}')
print('### time step 1')
stepper.advance()
t.assign(float(t) + float(dt))
U_restart = U.copy(deepcopy=True)
stage_restart = stepper.stages.copy(deepcopy=True)
print('### time step 2 (warm-up)')
stepper.advance()
U.assign(U_restart)
stepper.stages.assign(stage_restart)
stepper.solver.snes.setMonitor(snes_monitor)
print('### time step 2 (monitored)')
stepper.advance()
print('### RESULT cumulative KSP iterations after each Newton step: ' + ' '.join(str(c[2]) for c in cum)
      + ' | SNES norms: ' + ' '.join(f'{c[1]:.6e}' for c in cum))
