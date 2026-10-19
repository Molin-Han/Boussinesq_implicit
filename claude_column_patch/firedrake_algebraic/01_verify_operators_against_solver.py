'''
01  Check the hand-written operators of fa_common.py against the objects of the REAL solver.

The solver script Nonlinear_Boussinesq_Irk_SC.py is executed UNMODIFIED (its source is
read from the repository and exec'd up to, but not including, the time loop, i.e. up to
the line "# Time stepping"), one time step is taken with stepper.advance(), and then the
PETSc objects are opened:

    outer KSP -> HDivSchurPC (python) -> assembled auxiliary operator P (4 fields)
              -> PCFIELDSPLIT -> fieldsplit_1 KSP -> (S, A11)
              -> PCMG -> level smoothers -> firedrake.ASMStarPC / ASMVankaPC -> PCASM

Checks
  C1  A11 used by fieldsplit_1 == block (fields 0,1,2) of the assembled auxiliary operator
  C2  auxiliary operator re-assembled at the state of rest, block (0,1,2)
        == hand-written field-1 form  fa.field1_form            (to rounding)
  C3  pressure row of the auxiliary operator: no velocity coupling, (p,p) = (dt/2) delta M_p
  C4  -(2/dt) (velocity row, pressure column)^T == hand-written divergence matrix B
  C5  how far the operator actually used (linearised about the initial state, b != 0)
      is from the rest-state operator
  C6  the Schur complement of b in the hand-written field-1 operator
        == weighted H(div) form (u_x,w_x) + m (u_z,w_z) + c (div u, div w)
  C7  patches: the index sets of the PCASM of every MG level == fa.patch_index_sets
  C8  an independently rebuilt mesh/space (fa.build_hierarchy) gives the same matrices
  C9  MG level operators == hand-written form on the level meshes (rest state) and the
      Galerkin property P^T A_f P = A_c of the MG interpolation

usage:  python 01_verify_operators_against_solver.py [--patch star|vanka] [--nx 4] [--nz 6]
                                                     [--degree 2] [--dt 50] [--K 1]
'''
import os
import sys
import types
import argparse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

parser = argparse.ArgumentParser()
parser.add_argument("--patch", default="vanka", choices=["star", "vanka"])
parser.add_argument("--nx", type=int, default=4)
parser.add_argument("--nz", type=int, default=6)
parser.add_argument("--degree", type=int, default=2)
parser.add_argument("--dt", type=float, default=50.0)
parser.add_argument("--K", type=float, default=1.0)
parser.add_argument("--length", type=float, default=3.0e5)
parser.add_argument("--height", type=float, default=4.0e3)
args, _ = parser.parse_known_args()

import fa_common as fa
from fa_common import *          # noqa
from firedrake.dmhooks import get_function_space

delta = fa.shift_value(args.K, args.dt, args.length, args.height)
workdir = os.path.join(HERE, "work", f"run01_{args.patch}_deg{args.degree}")
os.makedirs(workdir, exist_ok=True)
os.chdir(workdir)

# ------------------------------------------------------------------------------------
# run the real script up to the time loop
# ------------------------------------------------------------------------------------
path = os.path.join(fa.REPO, "Nonlinear_Boussinesq_Irk_SC.py")
src = open(path).read()
marker = "# Time stepping"
assert src.count(marker) == 1
src = src.split(marker)[0]
modname = "boussinesq_solver_script"
mod = types.ModuleType(modname)
mod.__file__ = path
sys.modules[modname] = mod
script_argv = ["--nx", str(args.nx), "--nz", str(args.nz), "--length", repr(args.length),
               "--height", repr(args.height), "--degree", str(args.degree),
               "--dt", repr(args.dt), "--tmax", repr(2 * args.dt), "--shift", repr(delta),
               "--no_rotation", "--patch", args.patch]
old_argv = sys.argv
sys.argv = [path] + script_argv
try:
    exec(compile(src, path, "exec"), mod.__dict__)
finally:
    sys.argv = old_argv

out = []


def report(line):
    PETSc.Sys.Print(line)
    out.append(line)


def reldiff(X, Y):
    D = (X - Y)
    nX = abs(X).max() if X.nnz else 0.0
    nD = abs(D).max() if D.nnz else 0.0
    return nD / nX if nX > 0 else nD


report(f"# solver script arguments: {' '.join(script_argv)}")
report(f"# delta = {delta:.12e}, c = dt/(2 delta) = {fa.c_factor(args.dt, delta):.6e}, "
       f"m = {fa.m_factor(args.dt):.6e}")

stepper = mod.stepper
W = mod.W
stepper.advance()          # first time step, as in the script (j == 0)

snes = stepper.solver.snes
outer_pc = snes.ksp.pc
hpc = outer_pc.getPythonContext()
report(f"outer PC python context: {type(hpc).__name__}; inner PC type: {hpc.pc.getType()}")
report(f"Newton iterations {snes.getIterationNumber()}, cumulative FGMRES iterations "
       f"{snes.getLinearSolveIterations()} (first time step, serial)")

P_used = fa.petsc_to_scipy(hpc.P.petscmat)
fis = fa.field_offsets(W)
is_u, is_uy, is_b, is_p = fis
is1 = np.concatenate([is_u, is_uy, is_b])
report(f"field sizes (u, u_y, b, p): {[len(i) for i in fis]}")

ksps = hpc.pc.getFieldSplitSubKSP()
ksp0, ksp1 = ksps
S1, A11 = ksp1.getOperators()
report(f"fieldsplit_1: ksp {ksp1.getType()}, pc {ksp1.pc.getType()}, "
       f"operator type {S1.getType()}, preconditioning matrix type {A11.getType()}")
A11_used = fa.petsc_to_scipy(A11)

# C1
c1 = reldiff(P_used[is1][:, is1].tocsr(), A11_used)
report(f"C1  max|A11(fieldsplit_1) - P[(0,1,2),(0,1,2)]| / max|A11| = {c1:.3e}")

# ------------------------------------------------------------------------------------
# re-assemble the auxiliary operator at the state of rest
# ------------------------------------------------------------------------------------
U_save = mod.U.copy(deepcopy=True)
st_save = stepper.stages.copy(deepcopy=True)
mod.U.assign(0)
stepper.stages.assign(0)
hpc.update(outer_pc)
P_rest = fa.petsc_to_scipy(hpc.P.petscmat)

W1 = mod.V_2D * mod.Vy * mod.Vb
a1, bcs1 = fa.field1_form(W1, args.dt, delta)
A_hand = fa.assemble_scipy(a1, bcs1)
A_rest = P_rest[is1][:, is1].tocsr()
report(f"C2  max|P_rest[(0,1,2),(0,1,2)] - A_hand| / max|A_hand| = {reldiff(A_hand, A_rest):.3e}   "
       f"(n = {A_hand.shape[0]}, nnz hand {A_hand.nnz}, nnz solver {A_rest.nnz})")
a1b, _ = fa.field1_form(W1, args.dt, delta, div_of_full_u=False)
A_hand_b = fa.assemble_scipy(a1b, bcs1)
report(f"C2b same with div(u_xz) instead of div(u_xz + u_y j): "
       f"max|diff| / max|A| = {reldiff(A_hand, A_hand_b):.3e}")

# C3 pressure row
Ppu = P_rest[is_p][:, is1]
Mp = fa.assemble_scipy(TestFunction(mod.Pressure) * TrialFunction(mod.Pressure) * dx)
Ppp = P_rest[is_p][:, is_p].tocsr()
report(f"C3  max|P_rest[p,(u,u_y,b)]| = {abs(Ppu).max() if Ppu.nnz else 0.0:.3e} "
       f"(relative to max|P[p,p]| = {abs(Ppp).max():.3e});  "
       f"max|P[p,p] - (dt/2) delta M_p| / max = {reldiff((args.dt / 2) * delta * Mp, Ppp):.3e}")

# C4 divergence
B_hand = fa.assemble_scipy(fa.divergence_form(mod.V_2D, mod.Pressure))
bc_nodes = fa.bc_dofs_velocity(mod.V_2D)
free = np.setdiff1d(np.arange(mod.V_2D.dim()), bc_nodes)
Pup = P_rest[is_u][:, is_p].tocsr()
B_solver = (-(2.0 / args.dt) * Pup.T).tocsr()
report(f"C4  max|B_hand - (-(2/dt) P[u,p]^T)| / max|B| on the free velocity DoFs = "
       f"{reldiff(B_hand[:, free].tocsr(), B_solver[:, free].tocsr()):.3e}   "
       f"(B is {B_hand.shape[0]} x {B_hand.shape[1]}, {len(bc_nodes)} velocity DoFs constrained)")

# C5 distance between the operator used and the rest-state operator
report(f"C5  max|A11_used - A_rest| / max|A_rest| = {reldiff(A_rest, A11_used):.3e}   "
       f"(operator used by the solver is linearised about the last Newton iterate, b != 0)")
D = (A11_used - A_rest).tocsr()
nA = sp.linalg.norm(A_rest)
report(f"C5b Frobenius: |A11_used - A_rest|_F / |A_rest|_F = {sp.linalg.norm(D) / nA:.3e}")

# C6 Schur complement of b
f1 = fa.field_offsets(W1)
j_u, j_uy, j_b = f1
Auu = A_hand[j_u][:, j_u].toarray()
Aub = A_hand[j_u][:, j_b].toarray()
Abu = A_hand[j_b][:, j_u].toarray()
Abb = A_hand[j_b][:, j_b].toarray()
S_elim = Auu - Aub @ np.linalg.solve(Abb, Abu)
m = fa.m_factor(args.dt)
c = fa.c_factor(args.dt, delta)
av, bcv = fa.velocity_form(mod.V_2D, m, c)
S_form = fa.assemble_scipy(av, bcv).toarray()
report(f"C6  max|S_elim - S_form| / max|S_form| = {abs(S_elim - S_form).max() / abs(S_form).max():.3e}   "
       f"(S_form = (u_x,w_x) + m (u_z,w_z) + c (div u, div w), m = {m:.6e}, c = {c:.6e})")
report(f"C6b max|A[u,u_y]|, |A[u_y,u]|, |A[u_y,b]|, |A[b,u_y]| = "
       f"{abs(A_hand[j_u][:, j_uy]).max() if A_hand[j_u][:, j_uy].nnz else 0.0:.1e}, "
       f"{abs(A_hand[j_uy][:, j_u]).max() if A_hand[j_uy][:, j_u].nnz else 0.0:.1e}, "
       f"{abs(A_hand[j_uy][:, j_b]).max() if A_hand[j_uy][:, j_b].nnz else 0.0:.1e}, "
       f"{abs(A_hand[j_b][:, j_uy]).max() if A_hand[j_b][:, j_uy].nnz else 0.0:.1e}  (u_y decouples)")
report(f"C6c symmetry of S_form: max|S - S^T| / max|S| = {abs(S_form - S_form.T).max() / abs(S_form).max():.3e}; "
       f"smallest eigenvalue {np.linalg.eigvalsh(0.5 * (S_form + S_form.T)).min():.3e}")

# restore the state and the operator
mod.U.assign(U_save)
stepper.stages.assign(st_save)
hpc.update(outer_pc)

# ------------------------------------------------------------------------------------
# MG levels: operators, patches, interpolation
# ------------------------------------------------------------------------------------
mgpc = ksp1.pc
nlev = mgpc.getMGLevels()
mgtypes = {0: "multiplicative", 1: "additive", 2: "full", 3: "kaskade"}
report(f"PCMG: {nlev} levels, type {mgtypes.get(int(mgpc.getMGType()), mgpc.getMGType())}")
level_mats = []
level_spaces = []
for lev in range(nlev):
    sm = mgpc.getMGSmoother(lev)
    Al, Pl = sm.getOperators()
    Vl = get_function_space(sm.getDM())
    level_spaces.append(Vl)
    Ml = fa.petsc_to_scipy(Pl)
    level_mats.append(Ml)
    ncol_l = Vl.mesh().unique().num_cells()
    line = (f"  level {lev}: {ncol_l} columns, {Vl.dim()} DoFs, smoother ksp {sm.getType()} "
            f"(max_it {sm.getTolerances()[3]}), pc {sm.pc.getType()}")
    if sm.pc.getType() == "python":
        ctx = sm.pc.getPythonContext()
        asmpc = ctx.asmpc
        subksps = asmpc.getASMSubKSP()
        sizes_petsc = [k.getOperators()[1].getSize()[0] for k in subksps]
        kind = "star" if isinstance(ctx, ASMStarPC) else "vanka"
        mine = fa.patch_index_sets(Vl, kind)
        sizes_mine = [len(p) for p in mine]
        # the index sets the solver's own PC object builds: the same call as in
        # ASMPatchPC.initialize, with the solver's options re-inserted in the database
        with stepper.solver.inserted_options(), fa._quiet_firedrake():
            theirs = ctx.get_patches(Vl)
        lg = Vl.dof_dset.lgmap
        theirs = [np.asarray(lg.applyIS(i).indices) for i in theirs]
        same_sets = (len(mine) == len(theirs)
                     and all(np.array_equal(np.sort(a), np.sort(b)) for a, b in zip(mine, theirs)))
        # patch matrices held by the PCASM == restriction of the level matrix to my index
        # sets, in the order in which Firedrake lists the DoFs of the patch
        errs = []
        for k, p in zip(subksps, mine):
            Mj = fa.petsc_to_scipy(k.getOperators()[1]).toarray()
            ref = Ml[p][:, p].toarray()
            errs.append(abs(Mj - ref).max() / abs(ref).max())
        line += (f"; python ctx {type(ctx).__name__}, PC type {asmpc.getType()} (asm.py sets PCASM type BASIC), {len(subksps)} patches, "
                 f"sub pc {subksps[0].pc.getType()}, patch sizes (PETSc) {sorted(set(sizes_petsc))}, "
                 f"(fa_common) {sorted(set(sizes_mine))}, identical index sets: {same_sets}, "
                 f"max|patch matrix of PCASM - level matrix restricted to my index set| / max = {max(errs):.2e}")
    report(line)

# C9 level operators against the hand-written form on the level meshes; Galerkin property
for lev in range(nlev):
    Vl = level_spaces[lev]
    # the level space of fieldsplit_1 is the mixed space (u, u_y, b) on the level mesh
    a_l, bcs_l = fa.field1_form(Vl, args.dt, delta)
    A_l = fa.assemble_scipy(a_l, bcs_l)
    report(f"C9  level {lev}: max|A_level(solver) - A_hand(rest)| / max|A_hand| = "
           f"{reldiff(A_l, level_mats[lev]):.3e}")
for lev in range(1, nlev):
    Pmat = mgpc.getMGInterpolation(lev)
    nf, nc = Pmat.getSize()
    x = PETSc.Vec().createSeq(nc)
    y = PETSc.Vec().createSeq(nf)
    Pd = np.zeros((nf, nc))
    for j in range(nc):
        x.set(0.0)
        x.setValue(j, 1.0)
        x.assemble()
        Pmat.mult(x, y)
        Pd[:, j] = y.array_r
    Af = level_mats[lev].toarray()
    Ac = level_mats[lev - 1].toarray()
    G = Pd.T @ Af @ Pd
    # compare away from the rows/columns of constrained DoFs
    Vc = level_spaces[lev - 1]
    fc = fa.field_offsets(Vc)
    bcn = fc[0][fa.bc_dofs_velocity(Vc.sub(0))]
    keep = np.setdiff1d(np.arange(nc), bcn)
    err = abs(G[np.ix_(keep, keep)] - Ac[np.ix_(keep, keep)]).max() / abs(Ac).max()
    report(f"C9b interpolation level {lev - 1} -> {lev}: matrix {nf} x {nc}; "
           f"Galerkin check max|P^T A_f P - A_c| / max|A_c| on unconstrained DoFs = {err:.3e} "
           f"(A_f, A_c: level operators used by the solver)")

# ------------------------------------------------------------------------------------
# C8 independent rebuild
# ------------------------------------------------------------------------------------
mh2 = fa.build_hierarchy(args.nx, args.nz, args.length, args.height, refinement=mod.args.refinement)
V2, Vy2, Vb2, Q2, W2, W12 = fa.build_spaces(mh2[-1], args.degree)
a2, bcs2 = fa.field1_form(W12, args.dt, delta)
A2 = fa.assemble_scipy(a2, bcs2)
B2 = fa.assemble_scipy(fa.divergence_form(V2, Q2))
report(f"C8  independent rebuild with fa.build_hierarchy: max|A - A_hand| / max|A| = {reldiff(A_hand, A2):.3e}, "
       f"max|B - B_hand| / max|B| = {reldiff(B_hand, B2):.3e}")
p_a = fa.patch_index_sets(W1, args.patch)
p_b = fa.patch_index_sets(W12, args.patch)
report(f"C8b patches of the rebuilt space identical to those of the solver space: "
       f"{all(np.array_equal(a, b) for a, b in zip(p_a, p_b)) and len(p_a) == len(p_b)}")

with open(os.path.join(HERE, f"results_01_verify_{args.patch}_deg{args.degree}.txt"), "w") as f:
    f.write("\n".join(out) + "\n")
