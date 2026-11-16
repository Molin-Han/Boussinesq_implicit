from firedrake import *
import numpy as np
from firedrake.output import VTKFile
from petsc4py import PETSc
print = PETSc.Sys.Print
import utils
import ufl
from argparse import ArgumentParser
# from argparse_formatter import DefaultsAndRawTextFormatter

dt = Constant(1.0e+1)
shift = Constant(1.0e+1)

def vector_3D(u, uy):
    return (
        u + uy * utils.j()
    )


class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, u, v):
        # appctx = self.get_appctx(pc)
        W = u.function_space()
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        w = vector_3D(wxz, wy)
        # velo = uxz
        # w = wxz
        p = - Constant(1.0) / shift * div(velo)
        Jp = lhs(utils.SLB_velocity(velo, p, b, w, dt))
        Jp += lhs(utils.SLB_buoyancy(velo, b, q, dt))
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

    un, uny, bn, pn = split(Un)
    unp1, unp1ym, bnp1, pnp1 = split(Unp1)
    w_xz, wy, q, phi = TestFunctions(W)

    f = Function(V_2D).project(as_vector([sin(2 * pi * x), 0.0, sin(2* pi * z / height)]))
    g = Function(DG_km1).interpolate(sin(2*pi*x) + sin(2*pi*z))

    # DiricheletBC
    bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
    bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
    bcs = [bc1, bc2]

    u = vector_3D(un, uny)
    w = vector_3D(w_xz, wy)

    eqn = utils.SLB_velocity(u, pn, bn, w, dt)
    eqn -= inner(w, f) * dx
    eqn += utils.SLB_buoyancy(u, bn, q, dt)
    eqn -= g * q * dx
    eqn += utils.SLB_pressure(u, phi)
    shift_eqn = eqn + shift * pn * phi * dx
    Jp = derivative(shift_eqn, Un)

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
        'ksp_rtol': 1e-1,
        'ksp_view': ':kspview.txt',
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

    # nprob = NonlinearVariationalProblem(eqn, Un, bcs=bcs, Jp=Jp)
    nprob = NonlinearVariationalProblem(shift_eqn, Un, bcs=bcs)
    nsolver = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=params_schur)

    nsolver.solve()
    name = 'SLB_slice'
    file_slb = VTKFile(name+'.pvd')
    un, uny, bn, pn = Un.subfunctions
    un.rename("in-plane-vel")
    uny.rename("y-vel")
    bn.rename("buoyancy")
    pn.rename("pressure")
    f.rename("f_RHS")
    g.rename("g_RHS")
    file_slb.write(un, uny, bn, pn, f, g)

if __name__ == "__main__":
    nx=40
    length=1.0
    height=1.0
    nz=20
    solve_LB(nx, length, nz, height, shift, dt)
