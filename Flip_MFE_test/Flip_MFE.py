from firedrake import *

class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, u, v):
        print("HDivSchurPC.form() is being called!")
        W = u.function_space()
        velo, b = split(u)
        w, q = split(v)
        k = as_vector([0., 1.])
        Jp = (inner(velo, w) + div(velo) * div(w))*dx
        Jp += q * b *dx
        Jp -= div(w) * b * dx
        # Jp += div(velo) * q * dx
        _, bcs = super().form(pc, u, v)
        return (Jp, bcs)

mesh = PeriodicSquareMesh(3, 3, 1, direction='x', quadrilateral=True)
x, y = SpatialCoordinate(mesh)
V_2D = FunctionSpace(mesh, 'RTCF', 1)
Pressure = FunctionSpace(mesh, 'DG', 0)
Vb = FunctionSpace(mesh, 'DG', 0)
W = V_2D * Vb * Pressure

U = Function(W)
u, b, p = split(U)
w, q, phi = TestFunctions(W)

f = Function(V_2D).project(as_vector([sin(2*pi*x), sin(2*pi*y)]))
g = Function(Vb).interpolate(sin(2*pi*x) + sin(2*pi*y))

eqn = inner(w, u) * dx - div(w) * p * dx - div(w) * b * dx - inner(w, f) * dx
eqn += q * b * dx - g * q * dx
eqn += phi * div(u) * dx + p * phi *dx

pc_params = {
    # 'mat_view':':matSchurAuxPC.txt',
    'mat_view':'draw:tikz:matSchurAuxPC.tex',
    'pc_type':'ksp',
    'ksp_ksp_type': 'preonly',
    'ksp_pc_type':'lu',
    'ksp_ksp_monitor': None,
    'ksp_pc_factor_mat_solver_type':'petsc',
}
params_schur = {
    'mat_type': 'aij',
    'ksp_type': 'gmres',
    'snes_type':'ksponly',
    'snes_monitor': None,
    'ksp_monitor_true_residual': None,
    'mat_view':'draw:tikz:full_mat.tex',
    'pc_type': 'fieldsplit',
    'pc_fieldsplit_type': 'schur',
    'pc_fieldsplit_schur_fact_type': 'full',
    'pc_fieldsplit_0_fields': '2',
    'pc_fieldsplit_1_fields': '0,1',
    'fieldsplit_0': {
        'ksp_type': 'preonly',
        'pc_type': 'lu',
    },
    'fieldsplit_1': {
        'ksp_type': 'preonly',
        'ksp_monitor': None,
        'mat_view':':field_1_mat.txt',
        # 'mat_view':'draw:tikz:field_1_mat.tex',
        'pc_type': 'python',
        'pc_python_type': __name__ + '.HDivSchurPC',
        'helmholtzschurpc': pc_params,
        },
}
nprob = NonlinearVariationalProblem(eqn, U, bcs=None)
nsolver = NonlinearVariationalSolver(nprob, solver_parameters=params_schur)
nsolver.solve()