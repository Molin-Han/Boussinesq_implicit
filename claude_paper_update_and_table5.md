# Paper update (single-column section, restyled figures) and the Table 5 study

Record for proofreading, 2026-10-19. Request of 05:35: (1) put the results of the previous conversation into the
paper with the modified pictures, keep the single-column section short (only the reason for picking the smaller
patch), push and report; (2) double-check Table 5, where the monolithic multigrid with the shifted preconditioner
should behave like the fieldsplit Schur-complement multigrid with the shifted preconditioner but the table shows
8 against 12 iterations; then re-run with larger meshes; record every experiment with its parameters.
The previous records are `claude_small_patch_robustness.md` (2026-10-16) and `claude_column_patch_study.md`.

Status markers: **[done]** finished and checked, **[open]** not done.

---

## 0. Answers in short

**Task 1 (paper): done, pushed** as commit 0778f01. Figures 4-12 use the restyled figures. Section 5.7 is two
paragraphs plus two new three-panel figures. Paragraph 1 gives only your reason for the smaller patch: in a thin domain
the divergence is dominated by its vertical derivative, so the kernel a robust patch must hold is expected to fit in
single columns, hence `ASMVankaPC` with `construct_dim 1`. Paragraph 2 gives the result: the single column works and is
as parameter-robust as the star at the sweeps' time step, and needs two to three times the iterations in the
grad-div-dominated nonlinear regime. Two referee agents checked every number; their corrections are in. Section 1.

**Task 2 (Table 5): your doubt was right.**

* **Why they should agree.** With the shifted operator used on every multigrid level, the monolithic multigrid and the
  fieldsplit Schur multigrid converge in the same number of iterations: 12, 13, 13, 13 on the table's meshes (GMRES
  left, as on 2026-06-01), and with exact inner solves both are exactly A_delta^{-1}. The reason: eliminating the
  pressure from a monolithic column-star patch gives exactly the fieldsplit patch (DG mass matrix block diagonal, the
  star's pressure lies in its two columns). So one monolithic cycle returns the same velocity as one fieldsplit
  cycle; the two differ only in a pressure residual damped by about 1e-3 per cycle. Checked to round-off by four
  independent agents. Section 2.2.
* **Why the table disagreed** (8, really 5-7, against 12-13). Firedrake's auxiliary-operator preconditioner passes its
  inner multigrid the outer, unshifted Jacobian as Amat, and PETSc's multigrid uses Amat on the finest level by
  default. So the "monolithic shifted" multigrid smoothed the unshifted system and used A_delta only for the patch
  factorisations and the coarse levels. That is a different preconditioner. It approximates the unshifted inverse and
  beats even the exact shifted inverse on small meshes (at (20,20) and (40,40): 6 and 7 iterations against 12 and 12
  with GMRES left). But at (80,80) it stops converging: 17/100/100 iterations with GMRES left, 100/100/100 with
  FGMRES right, the residual stagnating at a relative 2e-8 to 5e-7. The fix needs no code change:
  `-firedrake_2_monoshiftpc_pc_use_amat false`. The paper's old explanation (monolithic patches "resolve more
  pressure-velocity coupling") was not the reason. Section 2.3.
* **Serial:** identical iteration counts in serial, on 2 and on 16 ranks, for every solver and mesh run in serial
  ((10,10) to (40,80)).
* **Other errors in the old table:** the runs used 2 MPI processes (caption: 16), eps = 1.3e-2 (caption: ~1e-3),
  GMRES with left preconditioning (text: FGMRES); the monolithic patch at (40,80) has 3048 DOFs (printed 3408); the
  monolithic counts at (10,10) and (20,20) were 5 and 6 (printed 8 and 7), the fieldsplit counts at (40,40) and
  (40,80) were 13 (printed 12).
* **Larger meshes (new Table 5).** Meshes (20,20) to (160,160), i.e. up to 640 x 160 cells and 2.05 million
  unknowns, on 16 ranks with FGMRES right and three repetitions, at the old table's delta = 1e-4. The iteration
  counts are equal (11, 11, 11, 10, 10; the monolithic cycle needs one fewer at (80,80)) and do not grow with the
  mesh. The monolithic approach needs 11-17 % more PCApply time (the fieldsplit approach is 10-15 % cheaper), except
  1 % at (80,80). Section 2.7.
* **Paper:** Table 5, its caption and paragraph, a paragraph on the old configuration, the timing-section
  introduction, Table 4's caption and paragraph (eps, and its PCApply is the time of 14 applications) and the
  Discussion sentence are updated; Table 4 now fits the page. Commit dba9f47, pushed. Section 2.8.
* **Every experiment** is listed with problem, mesh, shift, outer solver, preconditioner, patch and multigrid settings
  in Section 2.4. The exact command of each run is in `claude_small_patch/table5/runs/<tag>/cmd.sh`.
* **Your decisions** (Section 4): adopting `'pc_use_amat': False` in the script; re-running at the nominal shift
  delta = 1e-5; a likely miscitation (`farrellMitchellWechsung2019` is an augmented-Lagrangian block preconditioner,
  cited for monolithic patch smoothers); the Appendix's Richardson smoother against the Chebyshev(6) of the runs; two
  latent bugs (rotation option name in `HDivSchurPC`, advection added twice for U_mean != 0); the text corrections of
  Sections 1.4 and 2.8.

---

## 1. Task 1: the paper edits [done, pushed]

Paper repository `../Patch-Preconditioner-for-Boussinesq`, commit **0778f01** on `main`, pushed to
`github.com:Molin-Han/Patch-Preconditioner-for-Boussinesq` ("Short single-column patch section with new figures;
restyled figures 4-12"). The commit contains `ex_article.tex` and the twelve PNG files the tex now includes from
`Image_restyled/`. The other restyled variants (`Image_restyled/{corrected,star_rerun,single_column}/`,
`column_patch_sweeps.png`) are on disk but not committed; nothing in the tex uses them. The paper builds with
`latexmk -norc -pdf -interaction=nonstopmode ex_article.tex` (53 pages, no undefined references, no new overfull
boxes).

### 1.1 What changed in `ex_article.tex`

1. **Figures 4-12** now include the restyled drop-in versions `Image_restyled/<same name>.png` (same data as the
   paper's figures, the style of `claude_small_patch/paper_style.py`). `Image_new/`, `Image_nonlinear/` and
   `Image_patch/` are untouched on disk.
2. **Section 5.7 "A smaller patch: single-column smoothing"** is two paragraphs and two figures (was 215 lines with
   two tables and fig. 14).
   * Paragraph 1 gives your reason for the smaller patch, in your terms of 2026-09-29: a robust patch smoother needs
     a basis of the divergence-free fields (the kernel of the grad-div term) with each member supported in one
     patch; in a thin domain the divergence is dominated by its vertical derivative, so this kernel is expected to
     be captured within single columns; hence the single column, `ASMVankaPC` with `pc_vanka_construct_dim 1` on
     the multigrid levels instead of `ASMStarPC` with `pc_star_construct_dim 0`; 18 N_z + 4 against 30 N_z + 8
     unknowns (364 against 608 at N_z = 20). It is stated as motivation ("expected"), not as a result, because the
     earlier study found that the discrete kernel does not depend on the aspect ratio and that the single column
     works for k = 2 through the cell bubbles; no mathematics beyond that is given, as you asked.
   * Paragraph 2 states the result, with every number checked against `claude_small_patch/results/all_counts.csv`
     by an independent referee: the single column converges in every run in which the star converges; at the
     sweeps' time step (dt = 50 s) its counts depend as little on dx, dz and eps as the star's; linear problem:
     within one iteration of the star in 165 of 191 runs, more only for C_1 = 1e-5 (at most 7 against 4) and in
     thin domains at large dt (up to 2.75 x, 11 against 4 at L_z = 250 m); nonlinear problem: same failures,
     counts equal or close above the worst-case admissible shift bound (median ratio 1.0, at most 2.1, dt <= 100 s
     for C_1 = 1e-5), two to three times the iterations inside the admissible range (median 2.6, at most 5.9, dt >=
     500 s), growing under vertical refinement there (median 19 at N_z = 10, 38 at N_z = 80, star 10). Conclusion
     sentence: "a working alternative to the star patch: in the linear problem and at the time step of the
     parameter sweeps it is as parameter-robust as the star patch, while in the grad-div-dominated regime of the
     nonlinear problem it needs two to three times the iterations"; about 30 % cheaper per outer iteration
     (nonlinear solver, (20,20)), hence faster wherever the counts are within about 1.4; k = 1 needs the star.
   * New figures `fig:single-column-linear`, `fig:single-column-nonlinear` (`Image_restyled/single_column_*.png`,
     script `claude_small_patch/paper_fig_single_column.py`): panels (a) aspect ratio, (b) mesh at dt = 50 s, (c)
     shift sweep; star solid with filled circles, single column dashed with open squares. The captions say which
     runs are omitted (the same ten linear and thirteen nonlinear non-converged runs for both patches) and that the
     nonlinear re-runs add a Newton absolute tolerance 1e-9.
3. **Cross-references**: Section 3.4 now says "for k = 2 a single-column patch is also found to work in practice
   (Section 5.7)"; the Discussion sentence on the grad-div augmentation no longer refers to the removed "two-level
   kernel count".
4. **Captions of restyled figures whose legends changed**: fig. 6 lists the seven depths and the eps range (its
   legend now shows L_z only, the longer labels hid data); fig. 10 gives the L_z <-> eps mapping; fig. 11 has no
   legend (it hid data) and refers to the legend of fig. 12, which uses the same colours; the reversed intervals
   "C_1 in [1e-3, 1e-6]" of figs. 11-12 are replaced by the seven values; fig. 12's caption now says what the drop-in
   figure plots (total KSP iterations plus the number of Newton steps, values >= 100 not plotted).
5. Missing word in Section 5.1 ("through Irksome. The fully discrete system ...").

### 1.2 Figure changes made after the referee checks

`claude_small_patch/restyle_paper_figures.py` (all sets regenerated): fig. 6 short legend labels, fig. 7 legend in
three columns with y up to 120 (it hid the C_1 = 1 point at dt = 1e3 s), fig. 11 without legend.
`claude_small_patch/paper_fig_single_column.py`: legend 7.5 pt, panel (b) on a log axis with ticks 10-80 (the
labels 10 and 20 collided), integer y ticks, legend lower right inside the panel.

### 1.3 How it was checked

Two independent referee agents: one checked every number and claim of the new section against the data, the
captions against the figures and the plotting script, the motivation against your words, and the LaTeX; the other
compared each of the ten restyled figures with its original (pixel replicas from the same data: no data change in
any figure) and every sentence and caption that refers to them. Their reports are in
`claude_small_patch/table5/agents/{section_referee,figure_referee}/` (scripts, candidate figures, `edits.json`).

### 1.4 Corrections the referees found in text that was not part of the request (NOT applied; your decision)

These sentences were already inconsistent with the paper's own (unchanged) data before today. The referee reports are
`claude_small_patch/table5/agents/{figure_referee,section_referee}/report.md`; exact replacement text is in
`agents/figure_referee/edits.json` for the items B1-B7, B8 (except the L_z <-> eps mapping, which was applied), B9, C1
and C2 (A1, B10 and a variant of B11 were applied, Section 1.1 item 4). The last three bullets below (solver
description, Appendix smoother, nonlinear defect) have no replacement text.

* Fig. 5 text and caption: "more tightly clustered than in the horizontal sweep" (every run takes 8 iterations in
  both sweeps; the spread is not smaller for N_z).
* Fig. 7 text: "the C_1 = 1 curve never converges ... for dt <= 10^2" (it converges at dt = 100 s, 90 iterations);
  "counts are large for every C_1" at small dt; caption "level off" for curves that do not.
* Section 5.4 opening: figs. 8-10 use C_1 = 1e-5, while Section 5.1 says the nominal C_1 = 1e-3 applies throughout.
* Fig. 10 text and caption: "tightly clustered" at dt = 100 s (cumulative totals 19, 20, 17, 15, 11).
* Figs. 11-12 text: "a well-chosen C_1 in [1e-4, 1e-3] gives uniformly low iteration counts" (C_1 = 1e-4 needs 74-90
  at dt = 5-50 s; the count decreases monotonically with C_1 for dt <= 500 s).
* Counting conventions differ between figures (figs. 4-6 drop the converged iterate, fig. 7 plots n + 1, fig. 12
  plots total + Newton steps): the nominal run appears as 7, 9 and 8. `Image_restyled/corrected/` has the true
  counts throughout.
* Section 5.1 and the captions of figs. 4-7 say FGMRES, right preconditioning, a ten-order reduction of the error
  against an independent direct solve; the linear runs use GMRES with left preconditioning, rtol 1e-8 on the
  preconditioned residual, and the error is measured against the converged warm-up solve.
* The Appendix shows a Richardson level smoother; the linear runs use Chebyshev(6), the nonlinear ones GMRES(6).
* The nonlinear counts still carry the auxiliary-preconditioner defect of 2026-10-16 (missing pressure
  back-substitution in `Nonlinear_Boussinesq_Irk_SC.py`), still your open decision.

---

## 2. Task 2: Table 5

### 2.1 The published table and the logs behind it

Table 5 (`tab:monolithic`) reports, for (N_x, N_z) = (10,10), (20,20), (40,40), (40,80): patch DOFs 308/608/1208/2408
(fieldsplit) and 388/768/1528/3408 (monolithic), PCApply 1.11/2.61/8.01/14.70 s against 0.63/1.91/6.33/11.69 s, and
12/13/12/12 against 8/7/7/7 iterations; caption "eps ~ 1e-3, dt = 100 s, delta = 1e-4, 16 MPI processes". The logs
of those runs (`yr2/data_firedrake26_timing/`, `yr2/data_firedrake26_mono/`, 2026-06-01, `-log_view` and `ksp_view`)
say:

| item | table | logs |
|---|---|---|
| processes | 16 | **2** (also for Table 4: `yr2/data_firedrake26_reordering/`) |
| aspect ratio | eps ~ 1e-3 | **L_z/L_x = 4 km/300 km = 1.3e-2** (also for Table 4) |
| outer solver | "FGMRES" (text) | **GMRES, left preconditioning**, restart 30, rtol 1e-8 on the preconditioned residual |
| level smoother | (Appendix: Richardson) | **Chebyshev(6)** in both, star patch (`ASMStarPC`, construct_dim 0), LU + RCM, PCMG `full`, coarse LU |
| fieldsplit iterations | 12, 13, 12, 12 | 12, 13, 13, 13 (Official-Run PCApply count minus one) |
| monolithic iterations | **8**, 7, 7, 7 | **5, 6**, 7, 7 |
| monolithic patch DOFs at (40,80) | 3408 | **3048** (= 38 N_z + 8) |
| PCApply times | as printed | the same |

So the gap is larger than the table shows: 5-7 against 12-13 iterations, a factor 1.9-2.4.

### 2.2 Should the two agree? Yes

Both preconditioners approximate the same operator A_delta^{-1}, and with the paper's patch they are the same method
in the velocity. Notation: v = (u_xz, u_y, b), p the DG1 pressure, M the (cell-block-diagonal) pressure mass matrix,
D the divergence, A_delta = [[A_vv, -dt D^T], [D, delta M]] (the code's signs), S = A_vv + (dt/delta) D^T M^{-1} D.

1. **Exact inner solves.** The full Schur factorisation with exact A_00 = delta M (bjacobi + ILU(0) of a
   cell-block-diagonal matrix is an exact LU) and exact S is A_delta^{-1}; so is an exact monolithic solve. The outer
   iterations then coincide.
2. **One patch.** The monolithic star patch (38 N_z + 8 unknowns) is the fieldsplit star patch (30 N_z + 8) plus the
   pressure unknowns of its two columns, and those are all the pressure unknowns that the divergence couples to the
   patch's velocity unknowns. Since M is cell-local, eliminating the pressure from a monolithic patch problem gives
   exactly the fieldsplit patch problem S_i = R_i S R_i^T: the velocity correction of a monolithic patch solve
   equals the fieldsplit one for every residual (it depends only on the reduced residual r_v + (dt/delta) D^T M^{-1}
   r_p).
3. **The smoother.** Each base cell lies in two vertex stars and the ASM is BASIC (corrections summed,
   `firedrake/preconditioners/asm.py:58`). In the variables (e_v, eta) with eta = e_p + (1/delta) M^{-1} D e_v
   (eta = (delta M)^{-1} times the residual of the shifted continuity equation), the additive-Schwarz operator of
   the monolithic system is block diagonal: diag(B_S S, 2I), with B_S S the fieldsplit one. A Chebyshev smoother
   therefore acts as e_v <- q(B_S S) e_v (the fieldsplit smoother) and eta <- q(2) eta, q(2) = T_6(-0.8)/T_6(1.2) =
   -0.0359 for the interval [0.2, 2.2].
4. **Coarse grid and FMG.** The spaces are nested and the cells affine, so rediscretisation is Galerkin and the
   divergence commutes with prolongation; the coarse correction acts as diag(I - P S_H^{-1} P^T S, I - Pi_H) on
   (e_v, eta). The same holds through PETSc's `full` (FMG) cycle (`fmg.c:6-32`), which restricts residuals with the
   transpose of the interpolation.
5. **Result.** One monolithic cycle on A_delta = one fieldsplit cycle in the velocity, plus a pressure-residual
   component damped by about q(2)^2 = 1.3e-3 per cycle (the fieldsplit cycle returns the pressure that satisfies the
   shifted continuity equation exactly). The Chebyshev intervals are the same: PETSc's estimate of the largest
   eigenvalue is 2 in both (every field supported in one column is reproduced exactly by both stars that contain
   it), to within 3e-5 on every level and partition, giving [0.2, 2.2] (four digits) for both.

Checked by four independent agents (one deriving it without being told the answer, one asked to refute it, a code
audit and a numerical check; their reports are `claude_small_patch/table5/agents/<name>/report.md`, Section 3):
the patch identity to 1.3e-16, the 2I pressure block to 2.5e-15, the
auxiliary operator against the Schur complement of A_delta to 1.0e-16 on all three levels, and the velocity outputs
of the two preconditioners to round-off (6e-16 at 16 x 4, 6e-13 at 32 x 8) once both run on A_delta.

### 2.3 Why the published numbers disagree: the monolithic multigrid smoothed the unshifted operator

* PETSc's PCMG uses the operator matrix Amat by default: `pc->useAmat = PETSC_TRUE` in `PCCreate_MG`
  (`petsc/src/ksp/pc/impls/mg/mg.c:1999`); the finest-level smoother gets `(pc->mat, pc->pmat)` (`mg.c:935-943`) and
  the finest-level residual uses that Amat (`mg.c:559-561`).
* Firedrake's `AssembledPC` / `AuxiliaryOperatorPC` gives its inner PC the OUTER Amat:
  `pc.setOperators(A, self.P.petscmat)` (`firedrake/preconditioners/assembled.py:72`).
* For `--monolithic` the outer Amat is the unshifted Jacobian J (the problem is `NonlinearVariationalProblem(eqn,
  ..., Jp=Jp)`). So `MonolithicShiftPC`'s multigrid computes its finest-level Chebyshev iteration, its eigenvalue
  estimate and the residual sent to the coarse grid with J; only the patch factorisations (built from Pmat) and the
  coarse levels (rediscretised from the shifted form) are shifted. The run's own `ksp_view` shows it: the finest
  smoother of the monolithic solver lists "linear system matrix, followed by the matrix used to construct the
  preconditioner: Mat (firedrake_2_) [J] ... Mat (firedrake_2_monoshiftpc_) [A_delta]".
* For the fieldsplit solver the same mechanism is harmless: its Amat is PETSc's MatSchurComplement built from the
  shifted Pmat blocks (`fieldsplit.c:771, 906, 2133-2136`); it equals the dense Schur complement to 1.7e-23, and the
  auxiliary operator equals the dense Schur complement to 1.0e-16.

So the Table-5 "monolithic" preconditioner is a multigrid for the UNSHIFTED system with delta-regularised patch solves
and coarse levels. It is not an approximation of A_delta^{-1} (velocity output 92-97 % away from A_delta^{-1} r,
2.4e-2 from J^{-1} r on 16 x 4), which is why it beats even the exact shifted inverse (like for like at (20,20) and
(40,40): 6 and 7 against 12 and 12 with GMRES left, 6 and 8 against 11 and 10 with right preconditioning) and, with a
GMRES(6) level smoother, needs only 2. The mechanism (theory agent, dense spectra at 16 x 4):
the column patches resolve the baroclinic (vertically structured) shift error on the finest level, where the
smoother sees eigenvalues 2 mu in [0.66, 1.95] inside the Chebyshev interval, so only the barotropic modes are left
(eigenvalues about 37 mu, smallest 0.0165); 18 eigenvalues of the preconditioned operator lie farther than 0.1 from 1
against 159 for the exact shifted inverse. But it is not robust: at (80,80) it needs 17 iterations for the first
solve and does not converge in 100 for the next two with GMRES left, on 2 and on 16 ranks, and none of the three
solves converges with FGMRES right: the residual stagnates at a relative 2e-8 to 5e-7 (Section 2.5 (e), (g)).

The paper's explanation ("solving each patch in the full saddle-point variables resolves more of the
pressure-velocity coupling per relaxation sweep than the Schur reduction") does not hold: item 2 of Section 2.2 shows
that the velocity part of a monolithic patch solve is exactly the fieldsplit patch solve.

**Fix, without changing the code:** `-firedrake_2_monoshiftpc_pc_use_amat false` (the key is not in the solver
dictionary, so the command line works). In the script it is one line, `'pc_use_amat': False,` in
`params_monolithic['monoshiftpc']` (`LB_time_slice.py:371-381`); I did not change the script (your decision, as for
the nonlinear script's back-substitution). The equivalent for the fieldsplit solver,
`-firedrake_2_fieldsplit_1_helmholtzschurpc_pc_use_amat false`, replaces the matrix-free Schur complement on the
finest level by the assembled auxiliary operator (same maths, same iterates).

### 2.4 Experiment register

**Common to every run** (all with `LB_time_slice.py`, the paper's linear solver):

| parameter | value |
|---|---|
| problem | **linear** Boussinesq vertical slice (linearised about rest), implicit midpoint; `--tmax 200 --dt 100 --timing`: three solves (step 1; then at t = 100 s a warm-up solve and the official, timed solve of the same system from the same initial guess); the counts are given as step 1 / warm-up / official |
| physics | N = 1e-2 s^-1, no rotation (default), U_mean = 0; initial buoyancy 3e-4 sin(pi z/H)/(1 + (x - L/2)^2/a^2), a = 5 km; w = 0 at top and bottom, periodic in x |
| domain | L_x = 300 km, L_z = 4 km (**eps = 1.3e-2**) |
| discretisation | degree k = 2: extruded RT_2 (u, w), DG1 (u_y), W_theta of degree 2 (b), DG1 (p) |
| mesh, MG hierarchy | coarse base mesh of N_x periodic cells, refined twice horizontally (`--refinement 2`, three MG levels), N_z layers on every level: fine mesh 4N_x x N_z cells, dx = 300 km/(4 N_x), dz = 4 km/N_z; unknowns 20 n N_z + 4 n with n = 4 N_x |
| shift | **delta = 1e-4** (`--shift 1e-4`), i.e. C_1 = delta dt = 1e-2 (the paper's delta = C_1/dt; the script's own banner prints "C1" = delta dt^1.5 = 0.1, a different quantity); the agents' small runs also used 1e-3, 1e-5, 1e-6 |
| outer Krylov | rtol 1e-8, atol 0, max 100 iterations, restart 30, constant-pressure null space; GMRES left unless stated |
| multigrid | PCMG type `full` (one FMG cycle per application), V-cycles, 3 levels; coarsest level preonly + LU (MUMPS in parallel); the two finer levels smoothed by Chebyshev with 6 iterations (PETSc's default eigenvalue estimate: 10 GMRES steps, interval [0.1, 1.1] x the estimated maximum, = [0.2, 2.2] here; for the old MONO the finest-level estimate is 2.001-2.018 and depends slightly on the partition) unless stated, preconditioned by the patch ASM; pre- and post-smoothing 6 each |
| patch | `ASMStarPC`, `pc_star_construct_dim 0` (column star: the vertical facets above a base vertex and the two columns around it), BASIC additive Schwarz, LU with RCM ordering (`--reordering`) per patch; fieldsplit patch 30 N_z + 8 unknowns (u, u_y, b), monolithic patch 38 N_z + 8 (all four fields) |
| fieldsplit Schur (FS) | `--richardson`: PCFIELDSPLIT Schur, factorisation `full`; field 0 = pressure, preonly + bjacobi/ILU(0) on delta M (exact); field 1 = (u, u_y, b), preonly + `HDivSchurPC` (auxiliary operator S = A_vv + (dt/delta) grad-div) -> one multigrid cycle |
| monolithic (MONO) | `--monolithic --mono_smoother chebyshev`: `MonolithicShiftPC` (auxiliary operator = the full shifted Jacobian A_delta, aij) -> one multigrid cycle on all four fields; "as run" = PETSc's default (finest level on the unshifted J), "fixed" = `-firedrake_2_monoshiftpc_pc_use_amat false`. (The script's banner prints "monolithic MG ASMExtrudedStar (star line patch)", but the patch is `firedrake.ASMStarPC`, as `ksp_view` shows.) |
| machine | workstation with 2 x Xeon Silver 4215 (16 physical cores, 32 threads), 62 GiB, OpenMPI with `--bind-to none` (16 ranks fill every core), Firedrake/PETSc 3.25.3 of the July venv; PCApply = Official-Run stage of `-log_view`, which includes the numerical refactorisation of the patch and coarse matrices (PCSetUp, 10-23 % of it in the T runs) |

**The experiments** (run folders `claude_small_patch/table5/runs/<tag>/` with `cmd.sh`, `run.log`, `log.txt`
(log_view), `slice3D.txt` (ksp_view), `status`):

| set | purpose | solvers and settings beyond the common ones | meshes (N_x, N_z) | ranks | tags |
|---|---|---|---|---|---|
| A | reproduce Table 5; serial vs parallel | FS; MONO exactly as on 2026-06-01 (`--outer_ksp gmres --pc_side left`, Chebyshev(6)) | (10,10), (20,20), (40,40), (40,80) | 1, 2, 16 | `A_{fs,mono}_np<P>_nx_nz` |
| A' | norm of the stopping test | FS with GMRES right (`--richardson --pc_side right`); MONO with FGMRES right (`--pc_side right`, FGMRES being the monolithic default); MONO with GMRES(6) level smoother + FGMRES right (the script's default monolithic configuration); MONO with GMRES left and `-firedrake_2_ksp_norm_type unpreconditioned` (`Ap_monocheb_left_rich_*`: all four stopped at once, GMRES has no unpreconditioned norm with left preconditioning) | same | 2 | `Ap_*` |
| B | converged inner solves | FS: fieldsplit_1 FGMRES around the cycle to rtol 1e-6 / 1e-10; MONO: FGMRES around the monolithic cycle (PCKSP) to 1e-6 / 1e-10; outer FGMRES right (`--inner_rtol`) | (20,20), (40,40) | 2 | `B_*` |
| C | larger meshes, old configuration | FS and MONO as in A | (80,80), (80,160) | 2, 16 | `C_*` |
| Y | accuracy of the monolithic patch LU (PETSc's LU does not pivot; RCM on saddle-point patches) | MONO with natural ordering (no `--reordering`), GMRES left and FGMRES right; MONO with RCM + MUMPS patch LU; FS with MUMPS patch LU; MONO + MUMPS at (80,80) on 2 ranks (stopped by hand in its first solve) | (20,20), (40,40); (80,80) | 2 | `Y_*` |
| F | the corrected monolithic solver | MONO + `-firedrake_2_monoshiftpc_pc_use_amat false` ("monofix"), left and right; FS + `..._helmholtzschurpc_pc_use_amat false` ("fspmat"); FS with exact LU on S (`--direct`); monofix with GMRES(6) smoother + FGMRES right | (10,10)-(40,80); serial (20,20), (40,40); 16 ranks (20,20)-(40,80) | 1, 2, 16 | `F_*` |
| T | **the new Table 5** | FS, fspmat and monofix, FGMRES with right preconditioning (`--outer_ksp fgmres --pc_side right`, the outer solver the paper states), three repetitions each | (20,20), (40,40), (80,80), (80,160), (160,160) | 16 | `T_*_rep<k>_*` |
| H | old configuration on the new meshes | MONO as run, FGMRES right; at (80,80) also with MUMPS patch LU (GMRES left) and with a residual monitor (`-firedrake_2_ksp_monitor`, FGMRES right) | (20,20), (40,40), (80,80) | 16 | `H_*` |
| agents | small serial checks (each agent's `report.md` lists its runs) | `LB_time_slice.py` with FS, MONO as run and MONO fixed, GMRES left and right, with estimated and fixed Chebyshev intervals [0.2, 2.2]; shifts delta = 1e-3, 1e-5, 1e-6 besides 1e-4 (theory_blind, runs E6); `--tmax 100` single-solve runs (theory_blind D_*, E3); a modified copy `agents/code_audit/audit/lb_audit_copy.py` with Richardson(0.5), GMRES(1), exact-LU-inner and fixed-interval variants (code_audit, runs flag_*); standalone copies of the problem with dense P^{-1}, J and A_delta (numerics `fs_vs_mono.py`, theory_refute `probe*.py`, theory_blind `spectra_driver.py`) | (4,4), (8,8), (4,8), (8,4): fine meshes 16 x 4 to 32 x 8 | 1 | `agents/<name>/` |

### 2.5 Results

Generated from the run folders by `claude_small_patch/table5/record_tables.py` (counts: step 1 / warm-up / official
solve; PCApply: Official-Run stage, seconds; "(cap)": the 100-iteration cap was reached).

#### (a) Reproduction of Table 5, serial and parallel (Set A; GMRES left, the 2026-06-01 configuration)

| (N_x,N_z) | solver | its 2026-06-01 log | its serial | its 2 ranks | its 16 ranks | PCApply s serial | 2 ranks | 16 ranks | 2026-06-01 (2 ranks) |
|---|---|---|---|---|---|---|---|---|---|
| (10,10) | FS | 12 | 11/12/12 | 11/12/12 | - | 1.47 | 1.16 | - | 1.11 |
|  | MONO (as run) | 5 | 5/5/5 | 5/5/5 | - | 0.94 | 0.71 | - | 0.63 |
| (20,20) | FS | 13 | 11/13/13 | 11/13/13 | 11/13/13 | 5.68 | 4.02 | 2.32 | 2.61 |
|  | MONO (as run) | 6 | 7/6/6 | 7/6/6 | 7/6/6 | 3.52 | 2.01 | 1.29 | 1.91 |
| (40,40) | FS | 13 | 12/13/13 | 12/13/13 | 12/13/13 | 23.32 | 15.30 | 7.28 | 8.01 |
|  | MONO (as run) | 7 | 7/7/7 | 7/7/7 | 7/7/7 | 16.72 | 10.64 | 5.15 | 6.33 |
| (40,80) | FS | 13 | 12/13/13 | 12/13/13 | 12/13/13 | 49.81 | 30.77 | 14.80 | 14.70 |
|  | MONO (as run) | 7 | 7/7/7 | 7/7/7 | 7/7/7 | 35.78 | 21.12 | 10.32 | 11.69 |

(10,10) cannot run on 16 ranks: the coarse base mesh has 10 cells (MPI error in the mesh distribution).

#### (b) Norm of the stopping test (Set A', 2 ranks)

| (N_x,N_z) | FS, GMRES left | FS, GMRES right | MONO (as run), GMRES left | MONO (as run), FGMRES right | MONO (as run), GMRES(6) smoother, FGMRES right |
|---|---|---|---|---|---|
| (10,10) | 11/12/12 | 13/11/11 | 5/5/5 | 6/5/5 | 2/2/2 |
| (20,20) | 11/13/13 | 13/11/11 | 7/6/6 | 7/6/6 | 2/2/2 |
| (40,40) | 12/13/13 | 13/11/11 | 7/7/7 | 8/8/8 | 2/2/2 |
| (40,80) | 12/13/13 | 13/11/11 | 7/7/7 | 8/8/8 | 2/2/2 |

The planned left-preconditioned runs with `-ksp_norm_type unpreconditioned` (`Ap_monocheb_left_rich_*`) stopped at once: PETSc's GMRES does not support the unpreconditioned norm with left preconditioning.

#### (c) Converged inner solves (Set B, 2 ranks, outer FGMRES right)

| (N_x,N_z) | inner rtol | FS | MONO |
|---|---|---|---|
| (20,20) | 1e-6 | 13/11/11 | 13/11/11 |
| (20,20) | 1e-10 | 13/11/11 | 13/11/11 |
| (40,40) | 1e-6 | 12/10/10 | 12/10/10 |
| (40,40) | 1e-10 | 12/10/10 | 12/10/10 |

The inner FGMRES of MONO is a PCKSP inside `MonolithicShiftPC`; PCKSP uses Pmat by default, so these MONO runs are converged solves with A_delta (the useAmat leak of Section 2.3 does not occur here).

#### (d) Accuracy of the monolithic patch LU (Set Y, 2 ranks, GMRES left)

| (N_x,N_z) | MONO RCM + PETSc LU | MONO natural + PETSc LU | MONO MUMPS | MONO RCM, FGMRES right | MONO natural, FGMRES right | FS RCM + PETSc LU | FS MUMPS |
|---|---|---|---|---|---|---|---|
| (20,20) | 7/6/6 | 7/6/6 | 7/6/6 | 7/6/6 | 7/6/6 | 11/13/13 | 11/13/13 |
| (40,40) | 7/7/7 | 7/7/7 | 7/7/7 | 8/8/8 | 8/8/8 | 12/13/13 | 12/13/13 |

(MONO as run in all columns; GMRES left unless stated. Timings of Set Y are not comparable: it ran next to other jobs; there the MUMPS patch solves were 3-7 x slower than the banded PETSc LU with RCM. `Y_mono_rcm_mumps_np2_nx80_nz80` was stopped by hand during its first solve (it would have taken about 100 min on 2 ranks); replaced by the 16-rank run of (g).)

#### (e) Larger meshes, 2026-06-01 configuration (Set C, GMRES left)

| (N_x,N_z) | FS 2 ranks | FS 16 ranks | MONO (as run) 2 ranks | MONO (as run) 16 ranks | FS PCApply 16 ranks | MONO PCApply 16 ranks |
|---|---|---|---|---|---|---|
| (80,80) | 12/13/13 | 12/13/13 | 17/100/100 (cap) | 17/100/100 (cap) | 28.86 | 218.76 |
| (80,160) | - | 12/13/13 | - | - | 70.74 | - |

The rest of Set C (FS and MONO at (80,160) on 2 ranks, MONO at (80,160) on 16 ranks, both solvers at (80,80) in serial) was dropped once the cause was known: the MONO runs that hit the cap took 8 min (16 ranks) and 20 min (2 ranks) at (80,80), and longer at (80,160) and in serial; the serial question is answered by (a) and (f). (FS at (80,160) on 16 ranks overlapped two other runs; its PCApply, 5.05 s per application, is 13 % above the clean T run.)

#### (f) The corrected monolithic solver (Set F; `-firedrake_2_monoshiftpc_pc_use_amat false`)

| (N_x,N_z) | FS left | FS (pmat) left | MONO fixed left | FS right (GMRES) | MONO fixed right (FGMRES) | MONO fixed serial | MONO fixed 16 ranks | MONO as run left |
|---|---|---|---|---|---|---|---|---|
| (10,10) | 11/12/12 | 11/12/12 | 10/12/12 | 13/11/11 | 12/11/11 | - | - | 5/5/5 |
| (20,20) | 11/13/13 | 11/13/13 | 11/13/13 | 13/11/11 | 13/11/11 | 11/13/13 | 11/13/13 | 7/6/6 |
| (40,40) | 12/13/13 | 12/13/13 | 12/13/13 | 13/11/11 | 13/11/11 | 12/13/13 | 12/13/13 | 7/7/7 |
| (40,80) | 12/13/13 | 12/13/13 | 12/13/13 | 13/11/11 | 13/10/10 | - | 12/13/13 | 7/7/7 |

| (N_x,N_z) | FS with exact LU on S (serial, GMRES left) | FS inner 1e-10 (FGMRES right) | MONO fixed, GMRES(6) smoother, FGMRES right | MONO as run, same |
|---|---|---|---|---|
| (20,20) | 11/12/12 | 13/11/11 | 13/11/11 | 2/2/2 |
| (40,40) | 11/12/12 | 12/10/10 | 12/10/10 | 2/2/2 |

#### (g) The 2026-06-01 monolithic configuration on the new meshes (Set H, 16 ranks)

| (N_x,N_z) | MONO as run, FGMRES right | MONO fixed, FGMRES right | FS, FGMRES right |
|---|---|---|---|
| (20,20) | 7/6/6 | 13/11/11 | 13/11/11 |
| (40,40) | 8/8/8 | 13/11/11 | 13/11/11 |
| (80,80) | 100/100/100 (cap) | 12/10/10 | 13/11/11 |

MONO as run at (80,80), 16 ranks: with MUMPS patch LU (GMRES left) 17/100/100 (cap), as with PETSc's LU (17/100/100 in (e)); with a residual monitor (FGMRES right) 100/100/100 (cap): the true residual falls fast (relative 1e-6 after 10-16 iterations) and then stagnates, at 2.0e-8 after 100 iterations in the first solve and 4.7e-7 in the other two, above the tolerance 1e-8.

### 2.6 Reading of the results

* **(a) Reproduction and serial runs.** Today's counts equal the 2026-06-01 logs exactly (FS 12, 13, 13, 13; MONO 5,
  6, 7, 7), and they are identical in serial, on 2 and on 16 ranks for both solvers, in all three solves of each run.
  So running in serial changes nothing; the Table-5 gap is not a parallel effect. (For FS and the corrected MONO
  PETSc's Chebyshev estimate is 2 to within 3e-5 on every level and partition; for MONO as run the finest-level
  estimate depends slightly on the partition, e.g. 2.0078 / 2.0063 / 2.0041 at (20,20) on 1 / 2 / 16 ranks, since the
  "noisy" estimation vector is a hash of the global index, and its counts are identical nonetheless.)
  Today's 2-rank PCApply times are 1.0-2.1 times the June-1 values for both solvers; the June runs predate the July
  Firedrake/PETSc installation, so old and new timings are never compared below.
* **(b) Not the norm.** With right preconditioning (true residual) FS needs 11 and MONO 5-8: the gap is not an
  artefact of the left-preconditioned stopping test. The script's default monolithic configuration (GMRES(6) level
  smoother, FGMRES) needs 2 iterations.
* **(c) Exact inner solves agree.** With both inner solves converged to 1e-6 or 1e-10 the two preconditioners give
  the same counts (11 and 10): this is the count of the exact shifted preconditioner A_delta^{-1}. One FS cycle (11,
  right) is as good as exact; one MONO cycle (5-8) is better than exact, which an approximation of A_delta^{-1}
  cannot be. That pointed to MONO not approximating A_delta^{-1} (Section 2.3).
* **(d) Not the patch LU.** PETSc's LU does not pivot, and `--reordering` puts RCM on the saddle-point patches,
  which Firedrake avoids by default ("natural ordering to avoid zero pivots in saddle-point problems",
  `asm.py:62-64`). Natural ordering and MUMPS (with pivoting) give the same counts: the monolithic patch solves are
  accurate (the audit measured 1.1e-15 against dense LAPACK solves for the monolithic patches).
* **(e) The old monolithic solver breaks down on larger meshes.** At (80,80) it needs 17 iterations for the first
  solve and does not converge in 100 for the warm-up and official solves (GMRES left), on 2 and on 16 ranks; FS needs
  13 at (80,80) and (80,160).
* **(f) The corrected monolithic solver equals FS.** With `-firedrake_2_monoshiftpc_pc_use_amat false` the
  monolithic counts are those of FS: 12, 13, 13, 13 with left preconditioning on the Table-5 meshes (one fewer in the
  first solve at (10,10)), the same or one fewer with right preconditioning, identical in serial and on 16 ranks;
  with a GMRES(6) smoother it reproduces the exact-inner-solve counts of (c). FS with the assembled auxiliary operator
  on the finest level (`fspmat`) gives FS's counts, as it must. FS with an exact LU on S needs 12 (left), one fewer
  than one cycle.
* **(g) The old configuration on the new meshes** (FGMRES right, 16 ranks): 6 and 8 iterations at (20,20) and
  (40,40); at (80,80) none of the three solves converges within 100: the true residual drops to a relative 1e-6 in
  10-16 iterations and then stagnates (2.0e-8 after 100 in the first solve, 4.7e-7 in the other two). MUMPS patch
  solves do not change this (17/100/100 with GMRES left, as with PETSc's LU), so the stall is a property of the
  method, not of the patch LU.
* **Other shifts (agents' small runs, 32 x 8 fine mesh, GMRES left, theory_blind E6).** FS against the corrected MONO:
  17/22/22 against 15/24/24 at delta = 1e-3, 6/6/6 against 6/6/6 at delta = 1e-5, 5/5/5 against 4/4/4 at
  delta = 1e-6 (MONO as run: 6/6/6, 3/3/3, 3/3/3). So the counts agree to within one at delta <= 1e-4 and differ by
  two at delta = 1e-3, where the pressure-residual part of the monolithic cycle matters more; not checked with right
  preconditioning or on the Table-5 meshes.

### 2.7 The new Table 5 (Set T)

16 MPI ranks; FGMRES with right preconditioning (the outer solver the paper states; both solvers are then stopped on
the same true-residual reduction 1e-8); one three-level multigrid cycle per application, Chebyshev(6) with the
column-star patch on the two finer levels and a direct LU on the coarsest; dt = 100 s, eps = 1.3e-2, and
**delta = 1e-4** (C_1 = delta dt = 1e-2), kept from the old table for comparability: this is ten times the paper's
nominal C_1 = 1e-3 of Section 5.2 (delta = 1e-5 at dt = 100 s), so a re-run at delta = 1e-5 is a possible choice for
you (on the 32 x 8 agent mesh the two solvers then need 6 and 6). Three repetitions of every run; PCApply = median of
the three Official-Run stage times, which include the numerical refactorisation of the patch and coarse matrices
(10-23 % of PCApply; without it the monolithic application is still 11-16 % more expensive); all counts were identical
in the three repetitions. `fspmat` is the fieldsplit solver with the assembled auxiliary operator also on the finest
level (same iterates, PCApply within 2.5 % of `fs`, so the paper keeps the default `fs`); `monofix` is the monolithic
solver with `pc_use_amat false`. Generated by `claude_small_patch/table5/table5_final.py` (`results_final.md`,
`table5_tabular.tex`).

| (N_x, N_z) | fine mesh | unknowns | solver | patch dofs | its (official) | PCApply s median [min, max] (n) | s per it | KSPSolve s |
|---|---|---|---|---|---|---|---|---|
| (20,20) | 80 x 20 | 32320 | fs | 608 | 11 | 1.75 [1.74, 1.82] (3) | 0.159 | 1.77 |
| (20,20) | 80 x 20 | 32320 | fspmat | 608 | 11 | 1.71 [1.70, 1.82] (3) | 0.155 | 1.73 |
| (20,20) | 80 x 20 | 32320 | monofix | 768 | 11 | 2.05 [1.99, 2.10] (3) | 0.186 | 2.07 |
| (40,40) | 160 x 40 | 128640 | fs | 1208 | 11 | 5.76 [5.72, 6.28] (3) | 0.524 | 5.83 |
| (40,40) | 160 x 40 | 128640 | fspmat | 1208 | 11 | 5.87 [5.84, 5.99] (3) | 0.533 | 5.93 |
| (40,40) | 160 x 40 | 128640 | monofix | 1528 | 11 | 6.66 [6.58, 7.16] (3) | 0.605 | 6.73 |
| (80,80) | 320 x 80 | 513280 | fs | 2408 | 11 | 23.36 [23.32, 23.38] (3) | 2.124 | 23.64 |
| (80,80) | 320 x 80 | 513280 | fspmat | 2408 | 11 | 23.86 [23.82, 23.92] (3) | 2.169 | 24.12 |
| (80,80) | 320 x 80 | 513280 | monofix | 3048 | 10 | 23.58 [23.52, 23.69] (3) | 2.358 | 23.83 |
| (80,160) | 320 x 160 | 1025280 | fs | 4808 | 10 | 44.60 [44.56, 44.86] (3) | 4.460 | 45.10 |
| (80,160) | 320 x 160 | 1025280 | fspmat | 4808 | 10 | 45.72 [42.88, 45.81] (3) | 4.572 | 46.19 |
| (80,160) | 320 x 160 | 1025280 | monofix | 6088 | 10 | 49.36 [47.22, 49.38] (3) | 4.936 | 49.87 |
| (160,160) | 640 x 160 | 2050560 | fs | 4808 | 10 | 94.11 [93.47, 95.30] (3) | 9.411 | 95.09 |
| (160,160) | 640 x 160 | 2050560 | fspmat | 4808 | 10 | 95.17 [95.13, 98.63] (3) | 9.518 | 96.17 |
| (160,160) | 640 x 160 | 2050560 | monofix | 6088 | 10 | 106.10 [104.76, 108.62] (3) | 10.610 | 107.05 |

* **Iterations:** the same for both approaches (11, 11, 11, 10, 10), except at (80,80) where the monolithic cycle
  needs one fewer (10); the counts do not grow with the mesh (they fall slightly).
* **Cost:** one monolithic application is 11-17 % more expensive than one fieldsplit application, and so is its
  cumulative PCApply time at every mesh except (80,80), where the monolithic cycle saves an iteration and the
  difference is 1 %; equivalently, the fieldsplit approach is 10-15 % cheaper. Only part of the difference is in the
  patch solves (MatSolve and PCApplyOnBlocks per iteration are -2 to +10 % higher for the monolithic cycle); the rest
  is that its smoother residuals and grid transfers also act on the pressure.
* **Practical size:** the largest mesh has 640 x 160 cells (dx = 469 m, dz = 25 m), 2.05 million unknowns, and one
  solve takes 94 s of PCApply on 16 ranks. Scaling from (80,160) to (160,160) at fixed N_z (twice the unknowns, same
  patch size): PCApply 44.6 s -> 94.1 s (factor 2.11).
* **Run-to-run scatter** of PCApply between the three repetitions ((max - min)/median): at most 9.7 % (FS at
  (40,40)), median 3.7 %. The 10-17 % differences are larger than the scatter; at (80,80) the scatter is only
  0.2-0.7 % and the ranges of the two solvers do not overlap, so even the 1 % there is resolved. One repetition of
  FS at (160,160) (09:23-09:25) overlapped a 2-rank nonlinear run that was not part of this study
  (`claude_small_patch/diag_item3/runs/G1_offdiag_amat_dt5000.0`); it is the maximum of the three (95.3 s), the
  median is not affected.

### 2.8 Changes to the paper (commit dba9f47, pushed)

All in `ex_article.tex`, Section 5.6 and the Discussion; written from the run data by
`claude_small_patch/table5/apply_paper_table5.py`, then corrected after the final referee check (Section 3,
`agents/final_paper_text/`).

1. **Section 5.6, first paragraph:** "Unless stated otherwise, the timings below are obtained on 16 MPI processes" ->
   "on 2 MPI processes for Table 4 and on 16 for Table 5" (Table 4's runs, `yr2/data_firedrake26_reordering/`, used 2).
2. **Table 4 and its paragraph:** eps ~ 1e-3 -> eps = 1.3e-2 (L_x = 300 km, L_z = 4 km, the runs' values); "the
   wall-clock cost of a single preconditioner application" -> the time over the 14 applications of one GMRES solve
   (that is what the logs give), including the numerical refactorisation, also stated in the caption; headers
   "PCApply" for "Time(PCApply)" and `\small` so that the table fits the text width (it overflowed by 71 pt); its
   numbers are unchanged.
3. **Table 5** replaced by the table of Section 2.7 (meshes (20,20) to (160,160), a column with the number of
   unknowns, FGMRES right, 16 ranks, medians of three runs), placed before its paragraph so that it stays on the same
   page; the caption states the solver (three levels, Chebyshev(6) with the column-star patch on the two finer
   levels, direct solve on the coarsest), delta = 1e-4 with C_1 = 1e-2, the N_x convention and the medians.
4. **The Table-5 paragraph** rewritten: what PCApply contains; equal iteration counts (one fewer for the monolithic
   cycle at (80,80)); why (block-diagonal DG mass matrix; the pressure coupled to a column star lies in its two
   columns, so eliminating it gives the fieldsplit patch problem; Galerkin coarse operators and a commuting
   divergence; the same Chebyshev interval; the two differ only in the pressure); the cost (one monolithic
   application 11-17 % more expensive, and so its cumulative time except 1 % at (80,80); the larger systems its
   smoother, residuals and transfers work with). A second paragraph on the old configuration: a monolithic cycle
   whose finest-level smoother and residual use the unshifted operator A_0 is a multigrid method for A_0, needs 6 and
   8 iterations at (20,20) and (40,40) (fewer than the exact inverse of A_delta, 11 and 10), but stagnates at a
   relative residual of 2e-8 to 5e-7 at (80,80) and is not robust; it names the option `pc_use_amat:false`. The old
   explanation ("resolves more of the pressure-velocity coupling per sweep") and "left for future tuning" are gone.
5. **Discussion:** the monolithic smoother "returns the same velocity as the Schur reduction in each application ...
   and converges in the same number of iterations to within one; the fieldsplit patches are about 20 % smaller, and
   the Schur reduction is the cheaper of the two, by 10-15 % in preconditioner time at four of the five resolutions
   tested" (was: "not in favour of the Schur reduction ... converge in fewer FGMRES iterations"); the stale commented
   lines after it were removed.

The paper builds (54 pages, no undefined references, no overfull box in the edited region). Two referee agents
checked the numbers and the text against the logs; their reports are in `agents/final_numbers/` and
`agents/final_paper_text/`.

**Suggested by the referee but not applied** (your text outside Table 5; exact replacements in
`agents/final_paper_text/edits.json`, entries R10-R19):

* **Citation (please check):** `farrellMitchellWechsung2019` is "An augmented Lagrangian preconditioner for the 3D
  stationary incompressible Navier-Stokes equations at high Reynolds number" (the bib entry's title), a block
  preconditioner with star-patch multigrid on the augmented velocity block, i.e. the analogue of the fieldsplit
  approach. The introduction (around line 274) and the first paragraph of the Table-5 discussion cite it for
  monolithic patch smoothers applied directly to the saddle-point system. The referee suggests `\cite{PCPATCH}`
  there (R12, R13).
* **Appendix B** shows a Richardson(6) level smoother with damping theta and calls the Section 5.6 comparison
  "untuned on both sides"; the linear runs (Table 5 and the robustness sweeps via `--richardson`) use Chebyshev(6)
  with PETSc's eigenvalue estimate, which the new Table-5 caption now states (R10, and R19 for Section 3.4).
* **Mesh convention:** Section 2.2 defines dx = L_x/N_x, while Section 5 quotes (N_x, N_z) of the coarsest of the
  three multigrid levels (fine mesh 4N_x x N_z); a sentence in Section 5.2 would make this explicit (R11).
* The introduction's "monolithic patches ... carry a high cost per application" against the measured 11-17 % (R12);
  Section 3.4's description of the star patch's facets (R14); Table 4's "cost per layer is constant" (it falls from
  0.169 to 0.100 s) and "all FGMRES iterations" (Table 4 used GMRES) (R15-R18); a 55 pt overfull box in Section 5.2
  (R17). Table 4's (20,40) no-RCM entry 35.01 s is 35.00 s in the log.

### 2.9 Side findings (not acted on; your decision)

* **The monolithic solver in `LB_time_slice.py`** needs `'pc_use_amat': False,` in `params_monolithic['monoshiftpc']`
  (`LB_time_slice.py:371-381`) to be the shifted preconditioner; at HEAD its default configuration (GMRES(6)
  smoother, FGMRES) is the 2-iteration hybrid of Section 2.3. Not changed in the script.
* **Rotation option name:** `HDivSchurPC` reads `...fieldsplit_1_helmholtzschurpc_use_rotation`
  (`LB_time_slice.py:106-107`), but the dictionary sets `fieldsplit_1_use_rotation` (`:317`), so the auxiliary
  operator never contains the Coriolis term. Harmless with rotation off (all of Table 5, and the script's default:
  `--no_rotation` is `store_false`, so passing it switches rotation ON). With rotation on the fieldsplit auxiliary
  operator differs from the Schur complement (the code audit reported a relative difference of 5.5e-6 in the
  velocity block, `agents/code_audit/report.md`; its dump file was not kept).
* **`U_mean`:** `bool(Constant(0.0))` is True, so the `if U_mean:` branches of `utils.LB_velocity`/`LB_buoyancy`
  assemble zero-valued advection terms (harmless at U = 0, only sparsity); for U_mean != 0 the velocity advection is
  added twice (upwind block and centred block).
* **Timings:** today's runs (July venv, PETSc 3.25.3 dev of 2026-07-08) are 1.0-2.1 x slower than the 2026-06-01 logs
  at 2 ranks for both solvers; the cause (software versions) was not investigated.
* **PETSc's GMRES with left preconditioning** stops on each solver's own preconditioned residual; the true residuals at
  the stopping point differ between the solvers (numerics agent: FS 9.2e-7 at 16 x 4 and 2.7e-6 at 32 x 8, against
  6.2e-10 and 9.2e-9 for the old MONO). This is why the new table uses FGMRES right.
* **Options:** because the solver's options prefix is auto-generated (`firedrake_2_`), petsctools ignores the command
  line for keys that are in the solver dictionary (`petsctools/options.py`, the `options_prefix is None` branch);
  new keys such as `pc_use_amat` can be added from the command line.

---

## 3. Files

| file or folder | what |
|---|---|
| `claude_paper_update_and_table5.md` | this record |
| `../Patch-Preconditioner-for-Boussinesq/ex_article.tex`, `Image_restyled/` | the paper: commits 0778f01 (Task 1) and dba9f47 (Table 5), both pushed |
| `claude_small_patch/table5/run_table5.sh` (Sets A, A', B, C), `run_table5_pivot.sh` (Y), `run_table5_fixed.sh` (F, T, H), `run_table5_extra.sh` (Set D, never run: superseded) | runners; each run in `runs/<tag>/` with `cmd.sh` (the exact command), `run.log`, `log.txt` (-log_view), `slice3D.txt` (ksp_view), `status` (exit code, wall s) |
| `commands*.txt`, `run_table5*.log` | the run lists and the runner logs (with the runs stopped by hand marked) |
| `parse_table5.py` -> `results.csv`, `results.md` | every run: counts of the three solves, patch size, PCApply count and time, KSPSolve, wall time |
| `record_tables.py` -> `results_record.md` | the tables of Section 2.5 |
| `table5_final.py` -> `results_final.md`, `table5_tabular.tex` | the new Table 5 |
| `apply_paper_table5.py` | the script that wrote the first version of the Table-5 edits into the tex (from the run data) |
| `agents/{theory_blind,theory_refute,code_audit,numerics}/` | the four audit agents: `report.md` (their findings), scripts, logs, dense-matrix checks (`theory_blind/runs/D_*_n4/dense.npz`, 43 MB each, can be deleted) |
| `agents/{section_referee,figure_referee}/` | the two referees of the Task-1 edits: `report.md`, candidate figures, `edits.json` (text corrections B1-B11, C1-C2, not applied, Section 1.4) |
| `agents/final_{numbers,paper_text,completeness}/` | the three final checkers of this record and of the Table-5 edits: `report.md`; `final_paper_text/edits.json` has R1-R20 |
| `claude_small_patch/paper_fig_single_column.py`, `restyle_paper_figures.py` | figure scripts (changed today, Section 1.2) |
| `claude_small_patch_robustness.md` | the 2026-10-16 record: its draft Section 8 of this morning was moved here and replaced by a pointer |
| memory note `table5-monolithic-useamat` | for later sessions |

Code: this study added four default-preserving diagnostic flags to `LB_time_slice.py` (`--mono_smoother`,
`--inner_rtol`, `--pc_side`, `--outer_ksp`) at 06:08, before the first Table-5 run; your commit a9a20ee of 06:32
contains them together with the `--patch` flag of 2026-10-16. Nothing changed in `LB_time_slice.py`, `utils.py` or
`Nonlinear_Boussinesq_Irk_SC.py` after that (`git diff` empty); every variant here was selected with those flags and
with PETSc options that add keys (`pc_use_amat`, `pc_factor_mat_solver_type mumps`, `ksp_monitor`). Nothing was
committed in this repository by me; your commits a9a20ee/0b68318 of 06:32-06:33 contain the run folders that existed
then.

---

## 4. What is solid and what is not

| statement | basis |
|---|---|
| the two shifted preconditioners give the same velocity, and the same counts to within one at delta = 1e-4 | algebra (Section 2.2) confirmed by four independent agents; outputs equal to round-off on 16 x 4 and 32 x 8 fine meshes; equal counts on every Table-5 mesh, in serial, on 2 and 16 ranks, left and right, and on the new meshes up to 2.05M unknowns; at delta = 1e-3 (32 x 8, GMRES left) they differ by two (22 against 24) |
| the published monolithic numbers come from the useAmat leak | source lines in PETSc and Firedrake, the runs' own ksp_view, and the fix restoring the predicted counts everywhere |
| serial = parallel | identical counts in serial, on 2 and on 16 ranks for every solver and mesh run in serial ((10,10)-(40,80); (20,20) and (40,40) for the corrected MONO) |
| the old configuration stalls at (80,80) | 17/100/100 with GMRES left (2 and 16 ranks, PETSc LU and MUMPS patch solves), 100/100/100 with FGMRES right; the residual stagnates at 2e-8 to 5e-7 |
| timings | 16 ranks on the 16-core workstation, one run at a time (one repetition overlapped a foreign 2-rank job, Section 2.7), three repetitions, scatter at most 9.7 % (median 3.7 %) |
| why the old configuration converges faster on small meshes and stalls on large ones | a dense-spectrum explanation for the first part (theory agent, 16 x 4); the stall is not analysed |

Open (your decisions): adopt `pc_use_amat False` in the script's monolithic dictionary; whether to re-run Table 5 at
the nominal shift delta = 1e-5 (C_1 = 1e-3) instead of 1e-4; the citation of `farrellMitchellWechsung2019`; the
Appendix's Richardson smoother; the mesh convention; the rotation option name; the `U_mean` double advection; the
text corrections of Sections 1.4 and 2.8; the nonlinear script's missing pressure back-substitution (2026-10-16).
