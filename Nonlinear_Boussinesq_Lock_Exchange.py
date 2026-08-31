'''
Two dimensional lock exchange gravity current in a vertical slice, in the
*incompressible Boussinesq* setting.

This is the test case of Hiester, Piggott & Allison, "The impact of mesh adaptivity
on the gravity current front speed in a two-dimensional lock-exchange", Ocean
Modelling 38 (2011) 1-21, doi:10.1016/j.ocemod.2011.01.003, whose configuration
follows the DNS of Haertel, Meiburg & Necker, JFM 418 (2000) 189-212.

Equations solved (kinematic pressure p, buoyancy b, no background stratification
so N = 0, non rotating):

    du/dt + (u.grad)u = -grad(p) + b k + nu * lap(u)
    db/dt + (u.grad)b =  0
    div(u)            =  0

in a closed channel 0 <= x <= L = 0.8 m, 0 <= z <= D = 0.1 m.  The lock gate sits
at x_c = L/2: the fluid to the left is dense (b = -g'/2) and the fluid to the right
is light (b = +g'/2), with reduced gravity g' = g * drho/rho0 = 0.01 m/s^2
(drho/rho0 = 1e-3).  Removing the gate at t = 0 drives a dense current rightwards
along the bottom and a light current leftwards along the top.

Boundary conditions, as in the reference: no normal flow everywhere, *no slip* on
the bottom and *free slip* (rigid lid) on the top, so that a single run gives both
the no slip and the free slip front speed.  The end walls are free slip; they are
only felt once the fronts arrive, which is after the measuring window closes.  The
molecular viscosity is nu = 1e-6 m^2/s and there is no scalar diffusivity (kappa = 0,
Schmidt number Sc = infinity), giving a Grashof number

    Gr = (u_b * (D/2) / nu)^2 = 1.25e6,   u_b = sqrt(g' D / 2) = 0.0224 m/s.

The diagnostic is the non dimensional front speed (Froude number) Fr = u_g / u_b of
each current, obtained by linear regression of the front position over the window in
which the fronts travel from 0.2 m to 0.3 m from the lock.  The DNS values of
Haertel et al. (2000) are Fr = 0.675 (free slip, top) and Fr = 0.574 (no slip, bottom).

The spatial discretisation, the implicit midpoint time discretisation via Irksome and
the shifted Schur complement preconditioner (IRKAuxiliaryOperatorPC) are exactly those
of Nonlinear_Boussinesq_Mountain.py; the only new ingredient is the symmetric interior
penalty form of the viscous term, which the mountain test case does not need.
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
    description='Nonlinear incompressible Boussinesq lock exchange (vertical slice).',
    formatter_class=ArgumentDefaultsHelpFormatter
)

# ! Parameter settings.  The default mesh is 256 x 64, i.e. the resolution of case Z64 of
# ! the structured fixed mesh simulations that this test case is usually compared against.
parser.add_argument('--nx', type=int, default=64, help='Number of columns of the coarsest mesh in the hierarchy.')
parser.add_argument('--nz', type=int, default=64, help='Number of layers to extrude (same on every level).')
parser.add_argument('--length', type=float, default=0.8, help='Horizontal length L of the channel (m).')
parser.add_argument('--height', type=float, default=0.1, help='Height D of the channel (m).')
parser.add_argument('--refinement', type=int, default=2, help='Levels of the multigrid.')
parser.add_argument('--degree', type=int, default=2, help='Order of the element.')
parser.add_argument('--dt', type=float, default=0.02, help='Time stepping parameter.')
parser.add_argument('--tmax', type=float, default=25.0,
                    help='Time period that we solve.  The fronts reach the end walls at '
                         'about t = 26 s (top) and t = 31 s (bottom), and the measuring '
                         'window closes at about t = 23 s, so 25 s is enough.')
parser.add_argument('--dumpt', type=float, default=0.5, help='Time between two vtk dumps.')
parser.add_argument('--shift', type=float, default=None,
                    help='Shift parameter delta for the shift preconditioner. delta has units of '
                         '1/(velocity*length) = s/m^2, so by default it is set to '
                         'shift_scale/(u_b*length), which reproduces delta=1e-4 on the unit box.')
parser.add_argument('--shift_scale', type=float, default=1.0e-4,
                    help='Dimensionless shift delta*u_b*length, only used when --shift is not given.')

# ! Lock exchange test case settings (Hiester, Piggott & Allison 2011, section 3).
parser.add_argument('--g_prime', type=float, default=1.0e-2,
                    help="Reduced gravity g' = g*drho/rho0 (m/s^2).  The buoyancy jump across "
                         "the lock gate is exactly g', i.e. b = -g'/2 on the dense side and "
                         "b = +g'/2 on the light side.")
parser.add_argument('--nu', type=float, default=1.0e-6,
                    help='Kinematic viscosity (m^2/s).  Set to 0 for an inviscid run, which '
                         'drops the interior penalty term and makes the bottom free slip.')
parser.add_argument('--interface_width', type=float, default=0.0,
                    help='Half width of a tanh smoothing of the initial buoyancy interface (m).  '
                         'The reference uses a discontinuous lock gate, which is the default (0).')
parser.add_argument('--ip_penalty', type=float, default=10.0,
                    help='Constant C in the symmetric interior penalty parameter '
                         'sigma = C * degree^2 / h of the viscous term.')

# ! Diagnostic settings.  The front position is the extent of the dense (resp. light) fluid
# ! in a thin strip along the bottom (resp. top) wall, and the Froude number is fitted over
# ! the window in which the fronts have travelled between fit_start and fit_end from the lock.
parser.add_argument('--front_depth_frac', type=float, default=0.1,
                    help='Thickness of the wall strip used to measure the front position, as a '
                         'fraction of the channel height D.')
parser.add_argument('--fit_start', type=float, default=0.2,
                    help='Distance travelled from the lock at which the Froude number fit starts (m).')
parser.add_argument('--fit_end', type=float, default=0.3,
                    help='Distance travelled from the lock at which the Froude number fit ends (m).')

# ! Test settings
parser.add_argument('--show_args', action='store_true', help='Print all the arguments when the script starts.')
parser.add_argument('--rtol', type=float, default=1.0e-6, help='Relative tolerance for the ksp of linear solver.')
parser.add_argument('--atol', type=float, default=1.0e-30,
                    help='Absolute tolerance for the ksp of linear solver.  The buoyancy jump is '
                         'g\' = 1e-2 m/s^2 and the velocities are O(1e-2) m/s, so the residuals '
                         'are tiny in absolute terms and the ksp has to be driven by the relative '
                         'tolerance.')
parser.add_argument('--maxit', type=int, default=150, help='Max iteration number for the first ksp of the linear solve.')
parser.add_argument('--timing', action='store_true', help='If true, run the code without monitoring and test for the time.')

args = parser.parse_known_args()
args = args[0]

if args.show_args:
    PETSc.Sys.Print(args)

use_rotation = False  # ! The lock exchange is non rotating.
monitor_run = not args.timing
solver_name = 'MG ASMStar'

nx = args.nx
nz = args.nz
length = args.length
height = args.height
deg = args.degree
ar = height / length
nx_fine = nx * 2 ** args.refinement
deltax = length / nx_fine
deltaz = height / nz

# ! Derived scales of the test case.  u_b is the buoyancy velocity and T_b the gravity
# ! current time scale, both built on the *half* depth D/2, following Haertel et al. (2000).
g_prime = args.g_prime
u_buoy = (g_prime * height / 2.0) ** 0.5
t_buoy = (height / (2.0 * g_prime)) ** 0.5
grashof = (u_buoy * (height / 2.0) / args.nu) ** 2 if args.nu > 0.0 else float('inf')
courant = u_buoy * args.dt / deltax
# ! delta is dimensional (s/m^2); the tuned unit box value 1e-4 corresponds to the
# ! dimensionless group delta*U*L = 1e-4, which is what we keep fixed here.  The lock
# ! exchange has no mean flow, so the buoyancy velocity plays the role of U.
shift_value = args.shift if args.shift is not None else args.shift_scale / (u_buoy * length)
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("Two dimensional lock exchange, Boussinesq (incompressible) version of")
print("Hiester, Piggott & Allison (2011), Ocean Modelling 38, 1-21,")
print("following the DNS of Haertel, Meiburg & Necker (2000), JFM 418, 189-212.")
print(f"Physical domain with length {length} and height {height}, aspect ratio {ar}.")
print(f"Number of elements in x direction {nx_fine} and z direction {nz}, "
      f"deltax {deltax}, deltaz {deltaz}.")
print(f"Time stepping parameters dt {args.dt}, total time {args.tmax} "
      f"= {args.tmax / t_buoy} gravity current time scales.")
print(f"Shift parameter delta {shift_value} s/m^2, augmented Lagrangian parameter 1/delta "
      f"{1.0/shift_value}, dimensionless delta*u_b*L {shift_value * u_buoy * length}.")
print(f"Lock gate at x_c = {length/2} m, buoyancy jump g' = {g_prime} m/s^2 "
      f"(b = {-g_prime/2} on the dense side, b = {g_prime/2} on the light side), "
      f"interface half width {args.interface_width} m.")
print(f"Buoyancy velocity u_b = sqrt(g'*D/2) = {u_buoy} m/s, "
      f"time scale T_b = sqrt(D/(2*g')) = {t_buoy} s.")
print(f"Kinematic viscosity nu = {args.nu} m^2/s, no scalar diffusivity, "
      f"Grashof number Gr = (u_b*(D/2)/nu)^2 = {grashof}.")
print(f"Advective Courant number u_b * dt / deltax = {courant}.")
print("No slip on the bottom, free slip rigid lid on the top, free slip end walls.")
print("Reference front speeds (Haertel et al. 2000 DNS): Fr = 0.675 free slip (top), "
      "Fr = 0.574 no slip (bottom).")
print(f"The code is running with {solver_name} solver for the Schur complement created.")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )


def viscous_form(u, w, mesh, nu, penalty, no_slip_measures=()):
    '''
    Symmetric interior penalty (SIPG) form of -nu * lap(u) for the velocity, which lives
    in RT x DG and is therefore only H(div) conforming: the normal component is continuous
    but the tangential one jumps, so the tangential jumps have to be penalised.

    No slip is imposed weakly (Nitsche) on the measures in no_slip_measures.  It has to be
    weak: DirichletBC on an extruded RT space constrains the *facet* dofs of the bottom,
    which are the normal (vertical) velocity only, so u.n = 0 is strong but the tangential
    velocity at the wall is untouched by the strong bc.  Boundaries with no term here get
    the natural condition, zero tangential stress, i.e. free slip.
    '''
    n = FacetNormal(mesh)
    # ! CellDiameter, not a precomputed Constant(deltax): the auxiliary operator of the
    # ! preconditioner is coarsened onto the multigrid hierarchy, where the cells are bigger,
    # ! and a frozen fine grid penalty would lose coercivity on the coarse levels.
    h = CellDiameter(mesh)
    sigma = penalty * Constant(deg ** 2)
    F = nu * inner(grad(w), grad(u)) * dx
    for dS_ in (dS_v, dS_h):
        F -= nu * inner(avg(grad(w)), outer(jump(u), n('+'))) * dS_
        F -= nu * inner(outer(jump(w), n('+')), avg(grad(u))) * dS_
        F += nu * sigma / avg(h) * inner(jump(w), jump(u)) * dS_
    for ds_ in no_slip_measures:
        F -= nu * inner(grad(w), outer(u, n)) * ds_
        F -= nu * inner(outer(w, n), grad(u)) * ds_
        F += nu * sigma / h * inner(w, u) * ds_
    return F


# ! Closed channel, so an IntervalMesh rather than the PeriodicIntervalMesh of the mountain
# ! test case.  The two ends of the interval keep the markers 1 (x = 0) and 2 (x = L) through
# ! the extrusion and the embedding, and they are the side walls of the slice.
distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = IntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=args.refinement)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (args.refinement+1), extrusion_type='uniform')
# ! No orography here, so the plain (x, z) -> (x, 0, z) embedding of utils is enough.  That
# ! map is linear, hence exact at any set of nodes, which is why it can use the default DG
# ! variant where the mountain test case has to ask for 'equispaced'.
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
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

# DirichletBC, zero normal flow through all four walls of the closed channel.  On the
# top and bottom this constrains the vertical velocity, on the end walls 1 and 2 the
# horizontal one.  The tangential velocity is free everywhere; the no slip condition on
# the bottom is added weakly by viscous_form.
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bc3 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), 1)
bc4 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), 2)
bcs = [bc1, bc2, bc3, bc4]

nu = Constant(args.nu)
viscous = args.nu > 0.0
# ! Only the bottom is no slip; the top is a free slip rigid lid and the end walls are free
# ! slip, so ds_b is the only Nitsche measure.
no_slip_measures = (ds_b,)

# Initial condition: the fluid is at rest, dense (b = -g'/2) to the left of the lock gate
# and light (b = +g'/2) to the right of it.
xc = Constant(length / 2)
b_jump = Constant(g_prime)
u0_slice, u0yic, b0ic, p0ic = U.subfunctions  # ! subfunction for data assignment

# ! Irksome (bc_type='DAE') imposes both the algebraic constraint div(u) = 0 and the
# ! Dirichlet data on the *stage values* only, so with the implicit midpoint rule
# ! (u^n + u^{n+1})/2 is constrained and hence div(u^{n+1}) = -div(u^n) and
# ! u^{n+1}.n|_bdy = -u^n.n|_bdy: an initial state that is not solenoidal, or that
# ! does not respect the no normal flow condition, is never cleaned up, it just
# ! flips sign every step.  Here the fluid starts at rest, so u0 = 0 satisfies both
# ! exactly and the Leray projection that the mountain test case needs is not required.
u0_slice.assign(0.0)
u0yic.assign(0.0)
p0ic.assign(0.0)
if args.interface_width > 0.0:
    b_expr = 0.5 * b_jump * tanh((x - xc) / Constant(args.interface_width))
else:
    b_expr = conditional(x < xc, -0.5 * b_jump, 0.5 * b_jump)
b0ic.interpolate(b_expr)
print(f"Initial divergence: ||div(u0)||_L2 = {norm(div(u0_slice))}, initial buoyancy range "
      f"[{COMM_WORLD.allreduce(b0ic.dat.data_ro.min(), op=MPI.MIN)}, "
      f"{COMM_WORLD.allreduce(b0ic.dat.data_ro.max(), op=MPI.MAX)}].")

# Equations
u = vector_3D(uxz, uy)
w = vector_3D(w_xz, wy)

eqn = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=use_rotation)
if viscous:
    eqn += viscous_form(u, w, mesh, nu, args.ip_penalty, no_slip_measures=no_slip_measures)
# ! N2 = 0: the lock exchange has no background stratification, all of the buoyancy is
# ! carried by the prognostic field b, and the buoyancy equation is pure advection.
eqn += utils.Nonlinear_buoyancy_Irk(b, q, u, n, N2=Constant(0.0))
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
        if viscous:
            F += viscous_form(u, w, mesh, nu, args.ip_penalty, no_slip_measures=no_slip_measures)
        F += utils.Nonlinear_buoyancy_Irk(b, q, u, n, N2=Constant(0.0))
        F += utils.Nonlinear_pressure_Irk(u, phi)
        F += delta * p * phi * dx

        #  Boundary conditions
        bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
        bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
        bc3 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), 1)
        bc4 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), 2)
        bcs = [bc1, bc2, bc3, bc4]
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
    'snes_type': 'newtonls',  # ! We have not added the EW trick.
    'ksp_type': 'fgmres',
    'ksp_pc_side': 'right',
    'snes_atol': args.atol,
    'snes_rtol': args.rtol,
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
    params_schur['ksp_view'] = ':Nonlinear_lock_exchange_slice3D.txt'
else:
    params_schur.update({
        'snes_monitor': None,
        'ksp_converged_rate': None,
        'ksp_monitor_true_residual': None,
    })


butcher_tableau = GaussLegendre(1)  # ! Implicit Midpoint Rule.

stepper = TimeStepper(eqn, butcher_tableau, t, dt, U, bcs=bcs,
                      solver_parameters=params_schur, appctx=appctx)

# Diagnostics: the front positions of the two currents, the energetics, and the cell wise
# advective Courant number.
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
Courant_num_form = dt * (2 * avg(unn * v_dg0) * (dS_v + dS_h) + unn * v_dg0 * (ds_tb + ds_v))
Courant_denom = assemble(One * v_dg0 * dx)
Courant = Function(DG0, name="Courant")

# ! Front position.  In a thin strip of thickness h_probe along a wall the current is a
# ! single tongue attached to the end wall it came from, so the *area* of dense fluid in
# ! the bottom strip divided by h_probe is exactly the x coordinate of the nose of the
# ! bottom (no slip) current, and likewise from the right for the top (free slip) one.
# ! This is a good deal more robust than chasing a contour once the Kelvin-Helmholtz
# ! billows start shedding, and it needs no point evaluation, so it is cheap in parallel.
# ! The strip is snapped to a whole number of layers so that it ends on a facet: an
# ! indicator that cuts through the middle of a cell is only integrated to quadrature
# ! accuracy, which biases the front position by a fraction of a cell.  Dividing by the
# ! *measured* area of the strip rather than by length*h_probe removes what is left.
n_probe_layers = max(1, int(round(args.front_depth_frac * nz)))
h_probe = n_probe_layers * deltaz
bottom_strip = conditional(lt(z, Constant(h_probe)), 1.0, 0.0)
top_strip = conditional(gt(z, Constant(height - h_probe)), 1.0, 0.0)
dense_area_form = bottom_strip * conditional(lt(bn, 0.0), 1.0, 0.0) * dx
light_area_form = top_strip * conditional(gt(bn, 0.0), 1.0, 0.0) * dx
bottom_strip_area = assemble(bottom_strip * dx)
top_strip_area = assemble(top_strip * dx)
print(f"Front position measured in a wall strip of {n_probe_layers} layers, "
      f"h_probe = {h_probe} m = {h_probe / height} D.")

# ! Boussinesq energetics.  The potential energy density is -b*z (the density anomaly is
# ! rho' = -rho0*b/g), so the release of the initial step in b is what feeds the kinetic
# ! energy; with kappa = 0 any drift of the total is numerical.
kinetic_form = 0.5 * inner(un, un) * dx
potential_form = -bn * z * dx


def global_max(f, op=MPI.MAX):
    local = f.dat.data_ro.max() if op is MPI.MAX else f.dat.data_ro.min()
    return f.comm.allreduce(local, op=op)


def front_positions():
    '''
    (x of the bottom, dense, rightward moving nose, x of the top, light, leftward one).
    '''
    x_bottom = length * assemble(dense_area_form) / bottom_strip_area
    x_top = length * (1.0 - assemble(light_area_form) / top_strip_area)
    return x_bottom, x_top


def froude_number(times, positions, x_start, direction):
    '''
    Least squares fit of the front speed over the window in which the front has travelled
    between fit_start and fit_end from the lock, non dimensionalised by the buoyancy
    velocity.  direction is +1 for the rightward (bottom) front and -1 for the leftward
    (top) one.  Returns (Fr, number of samples in the window).
    '''
    times = np.asarray(times)
    positions = np.asarray(positions)
    travelled = direction * (positions - x_start)
    window = (travelled >= args.fit_start) & (travelled <= args.fit_end)
    if window.sum() < 2:
        return float('nan'), int(window.sum())
    speed = np.polyfit(times[window], travelled[window], 1)[0]
    return speed / u_buoy, int(window.sum())


def update_diagnostics():
    w_diag.project(dot(un, utils.k()))
    Courant_num = assemble(Courant_num_form)
    Courant.dat.data[:] = Courant_num.dat.data_ro / Courant_denom.dat.data_ro


if monitor_run:
    name = 'Nonlinear_Boussinesq_lock_exchange'
    file_lb = VTKFile(name+'.pvd')
    update_diagnostics()
    file_lb.write(un, uny, bn, pn, w_diag, Courant)

# Time series of the front positions, for the Froude number fit and for replotting later.
t_series = [0.0]
x_bottom_0, x_top_0 = front_positions()
x_bottom_series = [x_bottom_0]
x_top_series = [x_top_0]
print(f"Initial front positions: bottom (dense) {x_bottom_0} m, top (light) {x_top_0} m, "
      f"lock gate at {length/2} m.")

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
    # ! The fronts are tracked every step, not every dump: the Froude number is a linear
    # ! regression over a window of only about a tenth of the run, so it wants all the
    # ! samples it can get, and the two area integrals are far cheaper than a solve.
    x_bottom, x_top = front_positions()
    t_series.append(float(t))
    x_bottom_series.append(x_bottom)
    x_top_series.append(x_top)
    if tdump > dumpt - 0.5*args.dt:
        if monitor_run:
            update_diagnostics()
            file_lb.write(un, uny, bn, pn, w_diag, Courant)
            # ! max|uy| is a health check: nothing in this test case forces the y direction,
            # ! so uy has to stay at round off.  A growing uy means the embedding of the
            # ! slice in 3D has gone wrong somewhere.
            print(f"Dump at t = {float(t)} ({float(t)/t_buoy} T_b), "
                  f"front x_bottom = {x_bottom} m ({x_bottom - length/2} m travelled), "
                  f"front x_top = {x_top} m ({length/2 - x_top} m travelled), "
                  f"max w = {global_max(w_diag)}, min w = {global_max(w_diag, op=MPI.MIN)}, "
                  f"max|uy| = {max(abs(global_max(uny)), abs(global_max(uny, op=MPI.MIN)))}, "
                  f"KE = {assemble(kinetic_form)}, PE = {assemble(potential_form)}, "
                  f"max Courant = {global_max(Courant)}.")
        tdump -= dumpt

print("Iterations", itcount, "its per step", itcount/stepcount)

# Front speeds.  The reference values are the DNS of Haertel et al. (2000): the top,
# free slip current travels at Fr = 0.675 and the bottom, no slip one at Fr = 0.574.
fr_bottom, n_bottom = froude_number(t_series, x_bottom_series, length/2, +1.0)
fr_top, n_top = froude_number(t_series, x_top_series, length/2, -1.0)
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print(f"Froude numbers fitted over the window {args.fit_start} m to {args.fit_end} m "
      "travelled from the lock:")
print(f"  no slip (bottom) front:   Fr = {fr_bottom} from {n_bottom} samples, "
      f"DNS 0.574, error {100.0 * (fr_bottom - 0.574) / 0.574} %.")
print(f"  free slip (top) front:    Fr = {fr_top} from {n_top} samples, "
      f"DNS 0.675, error {100.0 * (fr_top - 0.675) / 0.675} %.")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

if COMM_WORLD.rank == 0:
    np.savetxt('Nonlinear_Boussinesq_lock_exchange_fronts.txt',
               np.column_stack([t_series, x_bottom_series, x_top_series]),
               header=f'nx {nx_fine} nz {nz} dt {args.dt} nu {args.nu} '
                      f'u_b {u_buoy} T_b {t_buoy}\nt  x_front_bottom  x_front_top')
