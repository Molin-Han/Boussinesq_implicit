from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
# from argparse import ArgumentParser
# from argparse_formatter import DefaultsAndRawTextFormatter

dt = Constant(1.0)
shift = Constant(1.0)

def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )


class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, u, v):
        W = u.function_space()
        One = as_vector([1., 1., 1.])
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        unph = Constant(0.5) * (velo + One)
        bnph = Constant(0.5) * (b + Constant(1.))
        w = vector_3D(wxz, wy)
        pnp1 = - Constant(1.) / shift * div(velo)
        Jp = lhs(utils.LB_velocity(velo, One, unph, w, bnph, pnp1, dt))
        Jp += lhs(utils.LB_buoyancy(b, Constant(1.), q, unph, dt))
        #  Boundary conditions
        bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
        bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
        bcs = [bc1, bc2]
        return (Jp, bcs)


# TODO: add the test and monitor for plotting
def solve_LB(nx, length, nz, height, delta, dt, tmax=100.0, refinement_levels=2):
    distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 2)}
    m = PeriodicIntervalMesh(nx, length,distribution_parameters=distribution_parameters)
    mh = MeshHierarchy(m, refinement_levels=refinement_levels)
    hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * (refinement_levels+1), extrusion_type='uniform')
    new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
    mesh = new_mh[-1]
    finest_mesh_name = "finest"
    mesh.name = finest_mesh_name

    x, y, z = SpatialCoordinate(mesh)
    deg = 2
    V_2D = utils.extrude_RT(mesh, k=deg) # RT2
    DG_km1 = FunctionSpace(mesh, 'DG', deg-1) # DG1
    Vy = DG_km1
    Pressure = DG_km1
    Vb = utils.W_theta(mesh, k=deg)
    W = V_2D * Vy * Vb * Pressure

    Un = Function(W)
    Unp1 = Function(W)

    unxz, uny, bn, pn = split(Un)
    unp1xz, unp1y, bnp1, pnp1 = split(Unp1)
    w_xz, wy, q, phi = TestFunctions(W)

    # Initial condition
    xc = Constant(length/2)
    yc = Constant(length/2)
    a = Constant(5000)
    U = Constant(0.) # TODO: there is no background flow.
    # This is a 4 components function.
    u0_slice, u0yic, b0ic, p0ic = Un.subfunctions
    u1_slice, u1yic, b1ic, p1ic = Unp1.subfunctions
    b0ic.project(sin(pi*z/height)/(1+((x-xc)**2)/a**2))
    b1ic.project(sin(pi*z/height)/(1+((x-xc)**2)/a**2))
    print('===============================================')
    print('Initial condition has been interpolated')
    name = 'ic'
    file_lb = VTKFile(name+'.pvd')
    u0, u0y, b0, P0 = Un.subfunctions
    file_lb.write(u0, u0y, b0, P0)
    print("Save initial condition.")

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

    # TODO: should change the equation
    eqn = utils.LB_velocity(unp1, un, unph, w, bnph, pnp1, dt)
    eqn += utils.LB_buoyancy(bnp1, bn, q, unph, dt)
    eqn += utils.LB_pressure(unp1, phi)
    shift = eqn + delta * pn * phi * dx
    Jp = derivative(shift, Un)

    # Pressure Nullspace
    v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
    nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])

    # helmholtz_schur_pc_params = {
    #         'ksp_type': 'preonly',
    #         'ksp_max_its': 30,
    #         'pc_type': 'mg',
    #         'pc_mg_type': 'full',
    #         'pc_mg_cycle_type':'v',
    #         'mg_levels': {
    #             # 'ksp_type': 'gmres',
    #             'ksp_type':'richardson',
    #             # 'ksp_type': 'chebyshev',
    #             'ksp_richardson_scale': 0.2,
    #             'ksp_max_it': 2,
    #             # 'ksp_monitor':None,
    #             "pc_type": "python",
    #             "pc_python_type": "firedrake.ASMStarPC",
    #             "pc_star_construct_dim": 0,
    #             "pc_star_sub_sub_pc_type": "lu",
    #             # "pc_star_sub_sub_pc_type": "svd",
    #             # "pc_star_sub_sub_pc_svd_monitor": None,
    #         },
    #         'mg_coarse': {
    #             'ksp_type': 'preonly',
    #             'pc_type': 'lu',
    #         },
    #     }

    helmholtz_schur_pc_params = {
        'pc_type':'ksp',
        'ksp_ksp_type': 'preonly',
        'ksp_pc_type':'lu',
        # 'ksp_ksp_monitor': None,
    }

    params_schur = {
        # 'mat_type': 'aij',
        'ksp_type': 'gmres',
        'snes_type':'ksponly',
        'ksp_atol': 0,
        'ksp_rtol': 1e-6,
        'ksp_view': ':slice3D.txt',
        'snes_monitor': None,
        # 'ksp_monitor': None,
        'ksp_monitor_true_residual': None,
        'pc_type': 'fieldsplit',
        'pc_fieldsplit_type': 'schur',
        'pc_fieldsplit_schur_fact_type': 'full',
        # 'pc_fieldsplit_0_fields': '2',
        # 'pc_fieldsplit_1_fields': '0,1',
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
            'helmholtzschurpc': helmholtz_schur_pc_params,
            },
    }

    nprob = NonlinearVariationalProblem(eqn, Unp1, bcs=bcs, Jp=Jp)
    nsolver = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=params_schur)

    # TODO: should change the time stepping.
    # Time Stepping
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
    dumpt = dt.dat.data[0]
    tdump = 0.
    dtc.assign(dt)
    i = 0
    while t < tmax.dat.data[0] - 0.5 * dt.dat.data[0]:
        print(f"=======================================The solver is currently solving for time:{t}==========================")
        t += dt.dat.data[0]
        tdump += dt.dat.data[0]
        i += 1
        nsolver.solve()
        Un.assign(Unp1)
        if tdump > dumpt - dt.dat.data[0]*0.5:
            file_lb.write(un, uny, bn, pn)
            tdump -= dumpt


if __name__ == "__main__":
    nx = 200
    length=1.0e3
    height=1
    # length = 3.0e5
    # height = 1.0e2
    nz = 50
    tmax = Constant(15.)
    dtc = dt
    shift_c = shift
    solve_LB(nx, length, nz, height, shift_c, dtc, tmax=tmax)
