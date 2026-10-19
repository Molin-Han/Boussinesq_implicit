"""
Numerical check of the factorisation claim for the auxiliary-operator preconditioner
of Nonlinear_Boussinesq_Irk_SC.py (HDivSchurPC.getNewForm).

Claim: with the Irksome stage Jacobians (GaussLegendre(1), one stage, stage-derivative form)
    A_delta = J(eqn + delta*p*phi*dx)                       (shifted form, no substitution)
    P_aux   = J(F_hat), F_hat = same with p -> p - (1/delta) div(u) substituted EVERYWHERE
one has P_aux = A_delta T, T = I except the pressure-row / velocity-column block
    T_pu = -(1/delta) M_p^{-1} D,   D = matrix of phi*div(u)*dx,  M_p = pressure mass matrix.
Consequently P_aux^{-1} = T^{-1} A_delta^{-1}: same velocity (and buoyancy) as the shifted
solve, but pressure p_var = p + (1/delta) M_p^{-1} D u = ((dt/2) delta)^{-1} M_p^{-1} r_p.

Run serially:
  OMP_NUM_THREADS=1 /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/python check_aux_factorisation.py
Writes check_aux_factorisation.txt next to this file.  Does not touch any other file.
"""
import os
import sys
import time
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

REPO = "/home/mh2520/molin/phd/yr2/Boussinesq_implicit"
HERE = os.path.join(REPO, "claude_small_patch", "diag_item3")
OUT = os.path.join(HERE, "check_aux_factorisation.txt")
sys.path.insert(0, REPO)   # for utils.py (read only)

from firedrake import *          # noqa: E402,F403
from irksome import GaussLegendre, MeshConstant   # noqa: E402
from irksome.stage_derivative import getForm      # noqa: E402
import utils                                      # noqa: E402

_lines = []


def out(s=""):
    s = str(s)
    sys.stdout.write(s + "\n")
    sys.stdout.flush()
    _lines.append(s)


def flush():
    with open(OUT, "w") as f:
        f.write("\n".join(_lines) + "\n")


def fro(X):
    return np.linalg.norm(X, "fro") if X.ndim == 2 else np.linalg.norm(X)


def rel(X, Y):
    """||X - Y|| / ||Y|| (Frobenius); returns the absolute norm if ||Y|| == 0."""
    nY = fro(Y)
    if nY == 0.0:
        return fro(X - Y)
    return fro(X - Y) / nY


def dense(petscmat):
    B = petscmat.duplicate(copy=True)
    B.convert("dense")
    return np.array(B.getDenseArray(), dtype=float, copy=True)


# --------------------------------------------------------------------------------------
# Problem set-up: tiny version of the script (nx=4, nz=3, degree 2, no hierarchy)
# --------------------------------------------------------------------------------------
nx, nz = 4, 3
length, height = 3.0e5, 1.0e4      # Skamarock-Klemp-type domain used by the sweeps
deg = 2
amplitude = 3e-4

t_start = time.time()
out("check_aux_factorisation.py  --  numerical check of P_aux = A_delta T")
out(f"date: {time.strftime('%Y-%m-%d %H:%M:%S')}")
out(f"mesh: PeriodicIntervalMesh(nx={nx}, L={length}) extruded nz={nz}, H={height}, embedded in 3D; degree {deg}")

m = PeriodicIntervalMesh(nx, length)
try:
    mh = MeshHierarchy(m, refinement_levels=0)
    hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz], extrusion_type="uniform")
    mesh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)[-1]
    out("embedding: utils.high_dim_mesh_hierarchy on a one-level hierarchy (as in the script)")
except Exception as e:  # fall back to the explicit embedding
    ext = ExtrudedMesh(m, layers=nz, layer_height=height / nz)
    xe, ze = SpatialCoordinate(ext)
    coord_fs = VectorFunctionSpace(ext, "DG", 1, dim=3)
    new_coord = assemble(interpolate(as_vector([xe, 0, ze]), coord_fs))
    mesh = Mesh(new_coord)
    mesh.init_cell_orientations(utils.j())
    out(f"embedding: manual (hierarchy path failed with {e!r})")

x, y, z = SpatialCoordinate(mesh)
n = FacetNormal(mesh)
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, "DG", deg - 1)
Pressure = FunctionSpace(mesh, "DG", deg - 1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure
N = W.dim()
dims = [Wi.dim() for Wi in W]
out(f"W = V_2D * Vy * Vb * Pressure, dims {dims}, total {N}")

# field index sets (global dof numbering of the monolithic matrix)
idx = [np.asarray(iset.indices, dtype=int) for iset in W.dof_dset.field_ises]
names = ["u(V_2D)", "uy", "b", "p"]
for nm, ii in zip(names, idx):
    out(f"  field {nm:8s}: {len(ii)} dofs, global indices {ii.min()}..{ii.max()}"
        f" (contiguous: {bool(np.all(np.diff(ii) == 1))})")
iu, iuy, ib, ip = idx
iv = np.concatenate([iu, iuy, ib])     # all non-pressure dofs ("velocity-ish" rows)
assert sorted(np.concatenate(idx).tolist()) == list(range(N))

# Nonzero smooth state (script's initial buoyancy + a smooth velocity so that the
# advection Jacobian is non-trivial; p enters linearly so p0 is irrelevant).
U = Function(W)
u0s, uy0, b0, p0 = U.subfunctions
xc = Constant(length / 2)
a = Constant(5000)
b0.project(amplitude * sin(pi * z / height) / (1 + ((x - xc) ** 2) / a ** 2))
U0 = Constant(1.0)
u0s.project(as_vector([U0 * cos(pi * z / height) * (1 + 0.3 * sin(2 * pi * x / length)),
                       0.0,
                       0.2 * U0 * sin(pi * z / height) * cos(2 * pi * x / length)]))
uy0.project(0.1 * U0 * sin(pi * z / height) * sin(2 * pi * x / length))
p0.project(1.0e-3 * sin(2 * pi * x / length) * cos(pi * z / height))
out(f"state: ||u0||={norm(u0s):.3e} ||uy0||={norm(uy0):.3e} ||b0||={norm(b0):.3e} ||p0||={norm(p0):.3e}")

uxz, uy, b, p = split(U)
w_xz, wy, q, phi = TestFunctions(W)
u = uxz + uy * utils.j()
w = w_xz + wy * utils.j()

bcs = [DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top"),
       DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")]

butch = GaussLegendre(1)
MC = MeshConstant(mesh)
t = MC.Constant(0.0)
dtc = MC.Constant(1.0)
delta_c = Constant(1.0)
stages = Function(W)


def build_forms():
    eqn = utils.Nonlinear_velocity_Irk(u, w, b, p, n, use_rotation=False)
    eqn += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
    eqn += utils.Nonlinear_pressure_Irk(u, phi)
    F_delta = eqn + delta_c * p * phi * dx          # shifted form, no substitution
    p_hat = p - Constant(1.) / delta_c * div(u)      # the script's substitution
    F_hat = utils.Nonlinear_velocity_Irk(u, w, b, p_hat, n, use_rotation=False)
    F_hat += utils.Nonlinear_buoyancy_Irk(b, q, u, n)
    F_hat += utils.Nonlinear_pressure_Irk(u, phi)
    F_hat += delta_c * p_hat * phi * dx              # p already substituted (as in getNewForm)
    return F_delta, F_hat


F_delta, F_hat = build_forms()


def stage_jacobian(F):
    Fnew, bcnew = getForm(F, butch, t, dtc, U, stages, bcs=bcs)
    J = derivative(Fnew, stages)
    A = assemble(J, bcs=bcnew, mat_type="aij").petscmat
    return dense(A), bcnew


# Pressure mass matrix and divergence matrix D (pressure test, trial in W; only the
# V_2D (and possibly Vy) columns are non-zero).  Both independent of dt, delta, stages.
p_tr, phi_t = TrialFunction(Pressure), TestFunction(Pressure)
Mp = dense(assemble(p_tr * phi_t * dx, mat_type="aij").petscmat)
Mp_inv = np.linalg.inv(Mp)
D_mixed = dense(assemble(derivative(phi * div(u) * dx, U), mat_type="aij").petscmat)
D = D_mixed[ip, :]                        # (n_p x N)
K_mixed = dense(assemble(derivative(div(w) * div(u) * dx, U), mat_type="aij").petscmat)
out()
out("--- D = matrix of phi*div(u)*dx (rows: pressure, columns: all of W) ---")
out(f"  ||D||_F = {fro(D):.6e}")
out(f"  rows of the mixed assembly outside the pressure block: {fro(np.delete(D_mixed, ip, axis=0)):.3e} (expect 0)")
out(f"  ||D[:, uy]|| = {fro(D[:, iuy]):.3e}, ||D[:, b]|| = {fro(D[:, ib]):.3e}, ||D[:, p]|| = {fro(D[:, ip]):.3e} (expect 0: div(uy j)=0 on the x-z manifold)")
out(f"  cond(M_p) = {np.linalg.cond(Mp):.3e}")
gd = rel(K_mixed[np.ix_(iu, iu)], D[:, iu].T @ Mp_inv @ D[:, iu])
out(f"  grad-div check ||K_uu - D^T M_p^-1 D|| / ||D^T M_p^-1 D|| = {gd:.3e}   (K = matrix of div(w)div(u)dx; 0 iff div(V_2D) is contained in Pressure)")
flush()


def analyse(delta, dt, stage_mode):
    out()
    out("=" * 100)
    out(f"CASE delta = {delta:g}, dt = {dt:g}, stage vector k = {stage_mode}")
    out("=" * 100)
    delta_c.assign(delta)
    dtc.assign(dt)
    if stage_mode == "zero":
        stages.assign(0.0)
    else:
        # smooth non-zero stage so that u0 + (dt/2) k differs from u0 (mid-Newton situation)
        stages.assign((0.5 / dt) * U)

    t0 = time.time()
    A, bcA = stage_jacobian(F_delta)
    P, bcP = stage_jacobian(F_hat)
    out(f"assembled both stage Jacobians ({N}x{N}) in {time.time() - t0:.1f} s")

    # Dirichlet dofs (stage bcs on W.sub(0) top/bottom)
    bc_nodes = np.unique(np.concatenate([np.asarray(bc.nodes, dtype=int) for bc in bcA]))
    bc_dofs = iu[bc_nodes]
    free = np.setdiff1d(np.arange(N), bc_dofs)
    # verify: identity rows and zeroed columns in both matrices (Firedrake's symmetric elimination)
    E = np.eye(N)
    row_err_A = fro(A[bc_dofs, :] - E[bc_dofs, :])
    row_err_P = fro(P[bc_dofs, :] - E[bc_dofs, :])
    col_err_A = fro(A[np.ix_(free, bc_dofs)])
    col_err_P = fro(P[np.ix_(free, bc_dofs)])
    out(f"Dirichlet dofs: {len(bc_dofs)} (of {len(iu)} V_2D dofs).  Identity-row error: A {row_err_A:.1e}, P {row_err_P:.1e}; "
        f"off-diagonal Dirichlet-column norm: A {col_err_A:.1e}, P {col_err_P:.1e}")
    out("  -> Firedrake zeroes the Dirichlet rows AND columns (diag 1), so the change of variables must not act on"
        " the constrained velocity dofs: T is built from D_bc = D with the Dirichlet columns zeroed (T_plain uses D as is).")

    D_bc = D.copy()
    D_bc[:, bc_dofs] = 0.0

    def make_T(Dmat):
        T = np.eye(N)
        T[ip, :] += -(1.0 / delta) * (Mp_inv @ Dmat)
        return T

    T_bc = make_T(D_bc)
    T_plain = make_T(D)
    Tinv_bc = 2.0 * np.eye(N) - T_bc          # T = I + L with L^2 = 0  =>  T^-1 = I - L
    out(f"  ||T_bc T_bc^-1 - I|| = {fro(T_bc @ Tinv_bc - np.eye(N)):.1e}")

    nA, nP = fro(A), fro(P)
    out(f"||A_delta||_F = {nA:.6e}, ||P_aux||_F = {nP:.6e}, cond(A_delta) = {np.linalg.cond(A):.3e}, cond(P_aux) = {np.linalg.cond(P):.3e}")

    # ---------------- matrix identity ----------------
    AT_bc = A @ T_bc
    AT_plain = A @ T_plain
    out()
    out("--- (1) matrix identity P_aux = A_delta T ---")
    out(f"  ||P_aux - A_delta T_bc||_F / ||A_delta||_F              = {fro(P - AT_bc) / nA:.3e}")
    out(f"  ||P_aux - A_delta T_bc||_F / ||P_aux||_F                = {fro(P - AT_bc) / nP:.3e}")
    out(f"  free dofs only: ||(P - A T_bc)[free,free]|| / ||A[free,free]|| = {fro((P - AT_bc)[np.ix_(free, free)]) / fro(A[np.ix_(free, free)]):.3e}")
    out(f"  with T_plain (D not restricted): ||P - A T_plain|| / ||A||     = {fro(P - AT_plain) / nA:.3e}"
        f"   [mismatch lives only in the Dirichlet columns: ||(P - A T_plain)[:, free]|| / ||A|| = {fro((P - AT_plain)[:, free]) / nA:.3e}]")

    # ---------------- block structure ----------------
    out()
    out("--- (2) block structure (Frobenius norms; blocks of the stage Jacobians) ---")
    Apu, Ppu = A[np.ix_(ip, iu)], P[np.ix_(ip, iu)]
    out(f"  pressure-row / V_2D-column block: ||P_pu|| = {fro(Ppu):.3e}, ||A_pu|| = {fro(Apu):.3e}, ratio = {fro(Ppu) / fro(Apu):.3e}"
        f"  (max|P_pu| = {np.abs(Ppu).max():.2e} vs max|A_pu| = {np.abs(Apu).max():.2e})")
    out(f"  pressure-row / uy-column block:   ||P_puy|| = {fro(P[np.ix_(ip, iuy)]):.3e}, ||A_puy|| = {fro(A[np.ix_(ip, iuy)]):.3e}")
    out(f"  pressure-row / b-column block:    ||P_pb|| = {fro(P[np.ix_(ip, ib)]):.3e}, ||A_pb|| = {fro(A[np.ix_(ip, ib)]):.3e}")
    App, Ppp = A[np.ix_(ip, ip)], P[np.ix_(ip, ip)]
    out(f"  A_pu vs (dt/2) D_bc[:,u]:                     rel diff = {rel(Apu, 0.5 * dt * D_bc[:, iu]):.3e}")
    out(f"  pressure/pressure: ||P_pp - A_pp|| / ||A_pp|| = {rel(Ppp, App):.3e};  ||P_pp - (dt/2) delta M_p|| / ||.|| = {rel(Ppp, 0.5 * dt * delta * Mp):.3e}")
    Aup, Pup = A[np.ix_(iu, ip)], P[np.ix_(iu, ip)]
    out(f"  velocity/pressure: ||P_up - A_up|| / ||A_up|| = {rel(Pup, Aup):.3e};  ||A_up - (-(dt/2) D_bc[:,u]^T)|| / ||A_up|| = {rel(Aup, -0.5 * dt * D_bc[:, iu].T):.3e}")
    out(f"  b/pressure and uy/pressure blocks: ||P_bp|| = {fro(P[np.ix_(ib, ip)]):.1e}, ||A_bp|| = {fro(A[np.ix_(ib, ip)]):.1e}, ||P_uyp|| = {fro(P[np.ix_(iuy, ip)]):.1e}, ||A_uyp|| = {fro(A[np.ix_(iuy, ip)]):.1e}")
    Auu, Puu = A[np.ix_(iu, iu)], P[np.ix_(iu, iu)]
    GD = (0.5 * dt / delta) * (D_bc[:, iu].T @ Mp_inv @ D_bc[:, iu])
    out(f"  velocity/velocity: ||P_uu - A_uu - (dt/2)(1/delta) D_bc^T M_p^-1 D_bc|| / ||(dt/2)(1/delta) D_bc^T M_p^-1 D_bc|| = {rel(Puu - Auu, GD):.3e}"
        f"   (||A_uu|| = {fro(Auu):.3e}, ||grad-div term|| = {fro(GD):.3e})")
    # every other block of P equals that of A
    others = [("uy", iuy), ("b", ib)]
    blocks = []
    for rn, ri in [("u", iu), ("uy", iuy), ("b", ib)]:
        for cn, ci in [("u", iu), ("uy", iuy), ("b", ib)]:
            if rn == "u" and cn == "u":
                continue
            blocks.append(f"{rn}{cn}:{rel(P[np.ix_(ri, ci)], A[np.ix_(ri, ci)]):.1e}")
    out("  remaining blocks, ||P_ij - A_ij|| / ||A_ij|| (abs if A_ij = 0): " + ", ".join(blocks))

    # ---------------- inverse / preconditioner action ----------------
    out()
    out("--- (3) preconditioner action: P_aux^-1 r vs T^-1 A_delta^-1 r for a random residual r ---")
    rng = np.random.default_rng(12345)
    r = rng.standard_normal(N)
    r[bc_dofs] = 0.0                      # stage residual vanishes on Dirichlet dofs
    xA = np.linalg.solve(A, r)
    xP = np.linalg.solve(P, r)
    out(f"  dense solve residuals: ||A xA - r||/||r|| = {fro(A @ xA - r) / fro(r):.1e}, ||P xP - r||/||r|| = {fro(P @ xP - r) / fro(r):.1e}")
    out(f"  velocity (V_2D): ||xP_u - xA_u|| / ||xA_u||   = {rel(xP[iu], xA[iu]):.3e}")
    out(f"  uy:              ||xP_uy - xA_uy|| / ||xA_uy|| = {rel(xP[iuy], xA[iuy]):.3e}")
    out(f"  buoyancy:        ||xP_b - xA_b|| / ||xA_b||   = {rel(xP[ib], xA[ib]):.3e}")
    p_formula = (Mp_inv @ r[ip]) / (0.5 * dt * delta)
    p_T = xA[ip] + (1.0 / delta) * (Mp_inv @ (D_bc @ xA))
    out(f"  pressure of P_aux^-1 r vs ((dt/2) delta)^-1 M_p^-1 r_p:          rel diff = {rel(xP[ip], p_formula):.3e}")
    out(f"  pressure of P_aux^-1 r vs p_A + (1/delta) M_p^-1 D_bc u_A:        rel diff = {rel(xP[ip], p_T):.3e}")
    out(f"  pressure of A_delta^-1 r vs ((dt/2) delta)^-1 M_p^-1 r_p:        rel diff = {rel(xA[ip], p_formula):.3e}   (NOT equal: shifted solve pressure depends on r_u)")
    out(f"  ||p_A|| = {fro(xA[ip]):.3e}, ||p_P|| = {fro(xP[ip]):.3e}, ||(1/delta) M_p^-1 D_bc u_A|| = {fro(p_T - xA[ip]):.3e}")
    # dependence on the velocity residual: change r_u only, keep r_p
    r2 = r.copy()
    r2[iv] = rng.standard_normal(len(iv))
    r2[bc_dofs] = 0.0
    xA2 = np.linalg.solve(A, r2)
    xP2 = np.linalg.solve(P, r2)
    out(f"  change r_u (keep r_p): pressure change of P_aux^-1 r: {rel(xP2[ip], xP[ip]):.3e};  of A_delta^-1 r: {rel(xA2[ip], xA[ip]):.3e}")
    # residual-based checks (do not rely on inverting P or A twice)
    x_pred = Tinv_bc @ xA
    out(f"  ||P (T_bc^-1 xA) - r|| / ||r|| = {fro(P @ x_pred - r) / fro(r):.3e}   (T^-1 A^-1 r solves P x = r)")
    out(f"  ||A (T_bc xP) - r|| / ||r||    = {fro(A @ (T_bc @ xP) - r) / fro(r):.3e}   (T P^-1 r solves A x = r)")
    Pinv = np.linalg.inv(P)
    Ainv = np.linalg.inv(A)
    out(f"  ||P_aux^-1 - T_bc^-1 A_delta^-1||_F / ||P_aux^-1||_F = {fro(Pinv - Tinv_bc @ Ainv) / fro(Pinv):.3e}")
    out(f"  ||P_aux^-1 - A_delta^-1||_F / ||A_delta^-1||_F         = {fro(Pinv - Ainv) / fro(Ainv):.3e}   (for scale: P_aux^-1 is NOT A_delta^-1)")
    # pressure rows of the inverse: P^-1[p, u-cols] should vanish (pressure output independent of r_u)
    out(f"  ||P_aux^-1[p, non-p cols]|| / ||P_aux^-1[p, :]|| = {fro(Pinv[np.ix_(ip, iv)]) / fro(Pinv[ip, :]):.3e};"
        f"  same for A_delta^-1: {fro(Ainv[np.ix_(ip, iv)]) / fro(Ainv[ip, :]):.3e}")
    out(f"  ||P_aux^-1[p, p] - ((dt/2) delta)^-1 M_p^-1|| / ||.|| = {rel(Pinv[np.ix_(ip, ip)], Mp_inv / (0.5 * dt * delta)):.3e}")
    flush()
    return dict(
        rel_identity=fro(P - AT_bc) / nA,
        rel_identity_free=fro((P - AT_bc)[np.ix_(free, free)]) / fro(A[np.ix_(free, free)]),
        pu_ratio=fro(Ppu) / fro(Apu),
        pp=rel(Ppp, App),
        up=rel(Pup, Aup),
        uu=rel(Puu - Auu, GD),
        p_formula=rel(xP[ip], p_formula),
        pA_formula=rel(xA[ip], p_formula),
        inv=fro(Pinv - Tinv_bc @ Ainv) / fro(Pinv),
    )


results = {}
results[(1e-6, 50.0, "zero")] = analyse(1e-6, 50.0, "zero")
results[(1e-9, 5000.0, "zero")] = analyse(1e-9, 5000.0, "zero")
results[(1e-6, 50.0, "nonzero")] = analyse(1e-6, 50.0, "nonzero")

out()
out("=" * 100)
out("SUMMARY (relative Frobenius norms)")
out("=" * 100)
out(f"{'delta':>8s} {'dt':>7s} {'k':>8s} | {'|P-AT|/|A|':>11s} {'free dofs':>10s} {'|P_pu|/|A_pu|':>13s} {'pp':>8s} {'up':>8s} {'uu grad-div':>11s} | {'p_P vs formula':>14s} {'p_A vs formula':>14s} {'|Pinv-Tinv Ainv|':>16s}")
for (d, dt, mode), R in results.items():
    out(f"{d:8.0e} {dt:7g} {mode:>8s} | {R['rel_identity']:11.2e} {R['rel_identity_free']:10.2e} {R['pu_ratio']:13.2e} {R['pp']:8.1e} {R['up']:8.1e} {R['uu']:11.2e} | {R['p_formula']:14.2e} {R['pA_formula']:14.2e} {R['inv']:16.2e}")
out()
out(f"total wall time {time.time() - t_start:.1f} s")
flush()
