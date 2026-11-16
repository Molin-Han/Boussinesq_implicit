from firedrake import *

class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, trial, test):
        print("HDivSchurPC.form() is being called!")
        W = trial.function_space()
        u, b = split(trial)
        w, q = split(test)
        Jp = inner(u, w)*dx
        Jp += q * b *dx
        Jp -= div(w) * b * dx
        _, bcs = super().form(pc, trial, test)
        return (Jp, bcs)

mesh = UnitSquareMesh(3, 3)
V_2D = FunctionSpace(mesh, 'RT', 1)
Vb = FunctionSpace(mesh, 'DG', 0)
W = V_2D * Vb

sol = Function(W)
U = TrialFunction(W)
u, b = split(U)
rhs_test = TestFunction(W)
w, q = split(rhs_test)
rhs = Function(W).assign(1.)

eqn = inner(w, u) * dx - inner(div(w), b) * dx + inner(q, b) * dx
L = inner(rhs, rhs_test) * dx

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
    'mat_view':'draw:tikz:field_1_mat.tex',
    'ksp_type': 'preonly',
    'ksp_monitor': None,
    'pc_type': 'python',
    'pc_python_type': __name__ + '.HDivSchurPC',
    'helmholtzschurpc': pc_params,
}
prob = LinearVariationalProblem(eqn, L, sol, bcs=None)
solver = LinearVariationalSolver(prob, solver_parameters=params_schur)
solver.solve()