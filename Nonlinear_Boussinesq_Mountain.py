'''
Nonhydrostatic flow past an "Agnesi" mountain in a vertical slice, in the
*incompressible Boussinesq* setting.

This is the Boussinesq analogue of the compressible Euler test case of section 3.3
of Cotter & Shipton, "A compatible finite element discretisation for the
nonhydrostatic vertical slice equations", GEM Int. J. Geomath. 14:25 (2023),
arXiv:2210.07861, whose reference implementation is
https://github.com/colinjcotter/sw_implicit/blob/master/slice_mountain_nh.py

Equations solved (kinematic pressure p, buoyancy perturbation b about the linear
background b_bg = N^2 z, f-plane rotation optional and off for this test case):

    du/dt + (u.grad)u + f x u + mu(z) (u.k) k = -grad(p) + b k
    db/dt + (u.grad)b + N^2 (u.k)             =  0
    div(u)                                    =  0

on a periodic slice -L/2 <= x <= L/2 (implemented as 0 <= x <= L), 0 <= z <= H,
with the lower boundary raised to the "witch of Agnesi" profile
z_s(x) = h_m a^2 / ((x - x_c)^2 + a^2) by the terrain following coordinate
transform z -> z + z_s (H - z)/H, and u.n = 0 on the top and bottom boundaries.
mu(z) is the Newtonian sponge of the paper, damping the vertical velocity above
z = z_B.  Initial data is a uniform wind u = (U, 0, 0) with b = 0.

The time discretisation is the implicit midpoint rule via Irksome, and the linear
systems are solved with the shifted Schur complement preconditioner
(IRKAuxiliaryOperatorPC) as in Nonlinear_Boussinesq_Irk_SC.py.
'''
from firedrake import *
from irksome import GaussLegendre, Dt, MeshConstant, TimeStepper, IRKAuxiliaryOperatorPC
import numpy as np
from mpi4py import MPI
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

parser = ArgumentParser(
    description='Nonlinear incompressible Boussinesq flow over an Agnesi mountain (vertical slice).',
    formatter_class=ArgumentDefaultsHelpFormatter
)

# ! Parameter settings
parser.add_argument('--nx', type=int, default=45, help='Number of columns of the coarsest mesh in the hierarchy.')
parser.add_argument('--nz', type=int, default=70, help='Number of layers to extrude (same on every level).')
parser.add_argument('--length', type=float, default=144.0e3, help='Horizontal length of our solution domain.')
parser.add_argument('--height', type=float, default=35.0e3, help='Height of our solution domain.')
parser.add_argument('--refinement', type=int, default=2, help='Levels of the multigrid.')
parser.add_argument('--degree', type=int, default=2, help='Order of the element.')
parser.add_argument('--dt', type=float, default=5.0, help='Time stepping parameter.')
parser.add_argument('--tmax', type=float, default=9000.0, help='Time period that we solve.')
parser.add_argument('--dumpt', type=float, default=500.0, help='Time between two vtk dumps.')
parser.add_argument('--shift', type=float, default=None,
                    help='Shift parameter delta for the shift preconditioner. delta has units of '
                         '1/(velocity*length) = s/m^2, so by default it is set to '
                         'shift_scale/(U_mean*length), which reproduces delta=1e-4 on the unit box.')
parser.add_argument('--shift_scale', type=float, default=1.0e-4,
                    help='Dimensionless shift delta*U_mean*length, only used when --shift is not given.')

# ! Mountain test case settings (Cotter & Shipton 2023, section 3.3, nonhydrostatic case)
parser.add_argument('--U_mean', type=float, default=10.0, help='Constant horizontal mean flow in x-direction (m/s).')
parser.add_argument('--half_width', type=float, default=1000.0, help='Half width a of the Agnesi mountain (m).')
parser.add_argument('--h_mount', type=float, default=1.0, help='Peak height of the Agnesi mountain (m).')
parser.add_argument('--mubar', type=float, default=0.15, help='Sponge strength, mubar = mu * dt (dimensionless).')
parser.add_argument('--z_sponge', type=float, default=25.0e3, help='Bottom z_B of the absorbing layer (m).')

# ! Test settings
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')
parser.add_argument('--rotation', action='store_true', help='If set, the Coriolis term is enabled (rotation ON, OFF for this test case).')
parser.add_argument('--rtol', type=float, default=1.0e-6, help='Relative tolerance for the ksp of linear solver.')
parser.add_argument('--atol', type=float, default=1.0e-30,
                    help='Absolute tolerance for the ksp of linear solver.  The mountain is a 1 m '
                         'perturbation of an exact steady state, so the residuals are tiny in '
                         'absolute terms and the ksp has to be driven by the relative tolerance.')
parser.add_argument('--maxit', type=int, default=150, help='Max iteration number for the first ksp of the linear solve.')
parser.add_argument('--timing', action='store_true', help='If true, run the code without monitoring and test for the time.')

args = parser.parse_known_args()
args = args[0]

if args.show_args:
    PETSc.Sys.Print(args)

use_rotation = args.rotation  # ! The nonhydrostatic mountain test case of the paper is non-rotating.
monitor_run = not args.timing
solver_name = 'MG ASMStar'

nx = args.nx
nz = args.nz
length = args.length
height = args.height
# ! delta is dimensional (s/m^2); the tuned unit box value 1e-4 corresponds to the
# ! dimensionless group delta*U*L = 1e-4, which is what we keep fixed here.
shift_value = args.shift if args.shift is not None else args.shift_scale / (args.U_mean * length)
deg = args.degree
ar = height / length
nx_fine = nx * 2 ** args.refinement
deltax = length / nx_fine
deltaz = height / nz
courant = args.U_mean * args.dt / deltax
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Flow over an Agnesi mountain, Boussinesq (incompressible) version of")
print("Cotter & Shipton (2023), section 3.3, nonhydrostatic regime.")
print(f"Physical domain with length {length} and height {height}, aspect ratio {ar}.")
print(f"Number of elements in x direction {nx_fine} and z direction {nz}, "
      f"deltax {deltax}, deltaz {deltaz}.")
print(f"Time stepping parameters dt {args.dt}, total time {args.tmax}.")
print(f"Shift parameter delta {shift_value} s/m^2, augmented Lagrangian parameter 1/delta "
      f"{1.0/shift_value}, dimensionless delta*U*L {shift_value * args.U_mean * length}.")
print(f"Mountain: half width a = {args.half_width} m, peak height h_m = {args.h_mount} m, "
      f"centred at x = {length/2} m.")
print(f"Mean flow U_mean = {args.U_mean} m/s in x-direction, "
      f"buoyancy frequency N = {float(sqrt(utils.buo_freq()))} 1/s.")
print(f"Nonlinearity parameter N*a/U = {float(sqrt(utils.buo_freq())) * args.half_width / args.U_mean} "
      "(=1 is the nonhydrostatic regime).")
print(f"Advective Courant number U_mean * dt / deltax = {courant}.")
print(f"Sponge layer above z_B = {args.z_sponge} m with mubar = mu*dt = {args.mubar}.")
print(f"Rotation is {'ON' if use_rotation else 'OFF'}.")
print(f"The code is running with {solver_name} solver for the Schur complement created.")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )


def mountain(x, xc, a, hm):
    '''
    "Mount Agnesi" (witch of Agnesi) orography profile z_s(x) = hm * a^2 / ((x-xc)^2 + a^2).
    '''
    return hm * a ** 2 / ((x - xc) ** 2 + a ** 2)


def mountain_mesh_hierarchy(mh, xc, a, hm, H, dim=3):
    '''
    Embed a hierarchy of 2D extruded slice meshes in dim dimensions and apply the
    terrain following coordinate transform (x, z) -> (x, z + z_s(x) * (H - z) / H)
    on every level of the hierarchy.  The map is analytic, so each level gets the
    same orography and the columns of the mesh are preserved (needed by the line /
    star smoother).  The levels are only nested up to the O(h^2) interpolation
    error of the transform, which is negligible here since hm / H ~ 3e-5.
    '''
    meshes = []
    for m in mh:
        x, z = SpatialCoordinate(m)
        zs = mountain(x, xc, a, hm)
        coord_fs = VectorFunctionSpace(m, "DG", 1, dim=dim)
        new_coord = assemble(interpolate(as_vector([x, 0, z + zs * (H - z) / H]), coord_fs))
        new_mesh = Mesh(new_coord)
        new_mesh.init_cell_orientations(utils.j())
        meshes.append(new_mesh)

    return HierarchyBase(meshes, mh.coarse_to_fine_cells,
                         mh.fine_to_coarse_cells,
                         mh.refinements_per_level, mh.nested)


distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=args.refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (args.refinement+1), extrusion_type='uniform')
xc = Constant(length/2)
a = Constant(args.half_width)
hm = Constant(args.h_mount)
new_mh = mountain_mesh_hierarchy(hierarchy, xc, a, hm, height, dim=3)
mesh = new_mh[-1]
finest_mesh_name = "finest"
mesh.name = finest_mesh_name

MC = MeshConstant(mesh)
dt = MC.Constant(args.dt)
shift = Constant(shift_value)
tmax = args.tmax
t = MC.Constant(0.0)
appctx = {
    "dt": dt,
    "shift": shift,
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

uxz, uy, b, p = split(U)  # ! split for writing the equations.
w_xz, wy, q, phi = TestFunctions(W)

# DiricheletBC, zero normal flow through the terrain following bottom and the flat top.
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]

# Sponge (Newtonian damping of the vertical velocity) in the top of the domain,
# mu(z) = 0 for z < z_B and mu(z) = mubar / dt * sin^2(pi/2 * (z-z_B)/(H-z_B)) above.
zB = Constant(args.z_sponge)
mubar = Constant(args.mubar)


def sponge(mesh, dtc):
    zm = SpatialCoordinate(mesh)[2]
    return conditional(zm <= zB, 0.0,
                       mubar / dtc * sin(0.5 * pi * (zm - zB) / (height - zB)) ** 2)


# Initial condition: a uniform horizontal wind U_mean, zero buoyancy perturbation
# (the background stratification N^2 z is carried by the linear term of the
# buoyancy equation) and zero pressure.
U_mean = Constant(args.U_mean)
u0_slice, u0yic, b0ic, p0ic = U.subfunctions  # ! subfunction for data assignment

# ! Irksome (bc_type='DAE') imposes both the algebraic constraint div(u) = 0 and the
# ! Dirichlet data on the *stage values* only, so with the implicit midpoint rule
# ! (u^n + u^{n+1})/2 is constrained and hence div(u^{n+1}) = -div(u^n) and
# ! u^{n+1}.n|_bdy = -u^n.n|_bdy: an initial state that is not solenoidal, or that
# ! does not respect the no normal flow condition, is never cleaned up, it just
# ! flips sign every step.  The initial data therefore has to satisfy both already.
# ! We therefore take the discrete Leray projection of the uniform flow onto the
# ! divergence free subspace of the RT space with u.n = 0 on top and bottom.
W_proj = V_2D * Pressure
U_proj = Function(W_proj)
u_pr, p_pr = TrialFunctions(W_proj)
w_pr, phi_pr = TestFunctions(W_proj)
a_proj = (inner(w_pr, u_pr) * dx - div(w_pr) * p_pr * dx + phi_pr * div(u_pr) * dx)
L_proj = inner(w_pr, as_vector([U_mean, 0., 0.])) * dx
bcs_proj = [DirichletBC(W_proj.sub(0), as_vector([0., 0., 0.]), "top"),
            DirichletBC(W_proj.sub(0), as_vector([0., 0., 0.]), "bottom")]
nullspace_proj = MixedVectorSpaceBasis(
    W_proj, [W_proj.sub(0), VectorSpaceBasis(constant=True, comm=COMM_WORLD)])
solve(a_proj == L_proj, U_proj, bcs=bcs_proj, nullspace=nullspace_proj,
      solver_parameters={'mat_type': 'aij',
                         'ksp_type': 'preonly',
                         'pc_type': 'lu',
                         'pc_factor_mat_solver_type': 'mumps',
                         'pc_factor_shift_type': 'inblocks'})
u0_slice.assign(U_proj.subfunctions[0])
print(f"Initial divergence: ||div(u0)||_L2 = {norm(div(u0_slice))}, "
      f"initial normal flow through the terrain: {assemble(abs(dot(u0_slice, n)) * ds_b)}.")

# Equations
u = vector_3D(uxz, uy)
w = vector_3D(w_xz, wy)

eqn = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=use_rotation)
eqn += sponge(mesh, dt) * inner(w, utils.k()) * inner(u, utils.k()) * dx
eqn += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
eqn += utils.Nonlinear_pressure_Irk(u, phi)

# Pressure Nullspace
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])


# Auxiliary Operator Preconditioner for Schur Complement Form of the Equation.
class HDivSchurPC(IRKAuxiliaryOperatorPC):
    _prefix = 'shiftedschurpc_'

    def getNewForm(self, pc, u0, test):
        prefix = (pc.getOptionsPrefix() or "") + self._prefix
        rotation = PETSc.Options().getBool(f"{prefix}use_rotation", False)  # ! Use this petsc option can make the equation to be consistent.
        appctx_PC = self.get_appctx(pc)  # ! returning dmhooks.get_appctx(pc).appctx.
        dtc = appctx_PC["dt"]
        delta = appctx_PC["shift"]
        n = appctx_PC["n"]
        W = u0.function_space()
        mesh = W.mesh()
        uxz, uy, b, p = split(u0)
        wxz, wy, q, phi = split(test)
        u = vector_3D(uxz, uy)
        w = vector_3D(wxz, wy)
        # ? The pressure elimination happened here, and no more pressure equation.
        p = p - Constant(1.) / delta * div(u)  # ! This gives the correct SC form and also the form needed for fieldsplit. Details in paper / notes.

        F = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=rotation)
        F += sponge(mesh, dtc) * inner(w, utils.k()) * inner(u, utils.k()) * dx
        F += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
        F += utils.Nonlinear_pressure_Irk(u, phi)
        F += delta * p * phi * dx

        #  Boundary conditions
        bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
        bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
        bcs = [bc1, bc2]
        return (F, bcs)


# ! The Correct IRKAuxOPPC Solver Parameter.
shifted_schur_pc_params = {
    'use_rotation': use_rotation,
    'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_fact_type': 'full',
    'pc_fieldsplit_schur_precondition': 'a11',
    'pc_fieldsplit_0_fields': '3',
    'pc_fieldsplit_1_fields': '0,1,2',
    'fieldsplit_0': {  # Doing a pure mass solve for the pressure block.
        'ksp_type': 'preonly',
        'ksp_reuse_preconditioner': None,  # ! double check if pc is reused in petsc log view. The pressure factorisation is the same all the time and will not need to be recomputed each time.
        'pc_type': 'python',
        'pc_python_type': 'firedrake.AssembledPC',
        'assembled_pc_type': 'bjacobi',
        'assembled_sub_pc_type': 'ilu',  # ILU needs an assembled matrix so that AssembledPC is needed.
    },
    'fieldsplit_1': {
        'ksp_type': 'preonly',  # ! need to tune this.
        # ! ASM line smoother on patch:
        'pc_type': 'mg',
        'pc_mg_type': 'full',
        'pc_mg_cycle_type': 'v',
        'mg_levels': {
            'ksp_type': 'gmres',
            'ksp_max_it': 6,  # ? more robust for larger max_it here.
            "pc_type": "python",
            "pc_python_type": "firedrake.ASMStarPC",
            "pc_star_construct_dim": 0,
            "pc_star_sub_sub_pc_type": "lu",
            'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
            'pc_star_sub_sub_pc_factor_reuse_ordering': None,
        },
        'mg_coarse': {
            'ksp_type': 'preonly',
            'pc_type': 'lu',
        },
    },
}


params_schur = {
    'mat_type': 'matfree',
    'snes_type': 'newtonls',
    'ksp_type': 'fgmres',
    'ksp_pc_side': 'right',
    'ksp_atol': args.atol,
    'ksp_rtol': args.rtol,
    'snes_max_it': 10,
    'ksp_max_it': args.maxit,
    'ksp_converged_maxits': None,
    'pc_type': 'python',
    'pc_python_type': __name__ + '.HDivSchurPC',
    'shiftedschurpc': shifted_schur_pc_params,
}

if args.timing:
    params_schur['ksp_view'] = ':Nonlinear_mountain_slice3D.txt'
else:
    params_schur.update({
        'snes_monitor': None,
        'ksp_converged_rate': None,
        'ksp_monitor_true_residual': None,
    })


butcher_tableau = GaussLegendre(1)  # ! Implicit Midpoint Rule.

stepper = TimeStepper(eqn, butcher_tableau, t, dt, U, bcs=bcs,
                      solver_parameters=params_schur, appctx=appctx)

# Diagnostics: vertical velocity (the field plotted in figure 4 of the paper) and
# the cell wise advective Courant number.
un, uny, bn, pn = U.subfunctions
un.rename("in-plane-vel")
uny.rename("y-vel")
bn.rename("buoyancy")
pn.rename("pressure")
w_diag = Function(Vb, name="w")

DG0 = FunctionSpace(mesh, "DG", 0)
One = Function(DG0).assign(1.0)
v_dg0 = TestFunction(DG0)
unn = 0.5 * (inner(-un, n) + abs(inner(-un, n)))  # gives fluxes *into* cell only
Courant_num_form = dt * (2 * avg(unn * v_dg0) * (dS_v + dS_h) + unn * v_dg0 * ds_tb)
Courant_denom = assemble(One * v_dg0 * dx)
Courant = Function(DG0, name="Courant")


def global_max(f, op=MPI.MAX):
    local = f.dat.data_ro.max() if op is MPI.MAX else f.dat.data_ro.min()
    return f.comm.allreduce(local, op=op)


def update_diagnostics():
    w_diag.project(dot(un, utils.k()))
    Courant_num = assemble(Courant_num_form)
    Courant.dat.data[:] = Courant_num.dat.data_ro / Courant_denom.dat.data_ro


if monitor_run:
    name = 'Nonlinear_Boussinesq_mountain'
    file_lb = VTKFile(name+'.pvd')
    update_diagnostics()
    file_lb.write(un, uny, bn, pn, w_diag, Courant)


# Time stepping
dumpt = args.dumpt
tdump = 0.
itcount = 0
stepcount = 0
while (float(t) < tmax - 0.5 * args.dt):
    tdump += args.dt
    if stepcount == 0:
        stepper.advance()
    else:
        with PETSc.Log.Stage('Solver'):
            stepper.advance()
    print(float(t))
    t.assign(float(t) + float(dt))
    itcount += stepper.solver.snes.getLinearSolveIterations()
    stepcount += 1
    if tdump > dumpt - 0.5*args.dt:
        if monitor_run:
            update_diagnostics()
            file_lb.write(un, uny, bn, pn, w_diag, Courant)
            print(f"Dump at t = {float(t)}, max w = {global_max(w_diag)}, "
                  f"min w = {global_max(w_diag, op=MPI.MIN)}, "
                  f"max Courant = {global_max(Courant)}.")
        tdump -= dumpt

print("Iterations", itcount, "its per step", itcount/stepcount)
