# Smaller (single-column) patch: the paper's robustness tests re-run, two ways to build the patch, and the dt trend of fig. 14(a)

Record for proofreading, written 2026-10-16 (session started 05:30). The 2026-10-19 follow-up is in `claude_paper_update_and_table5.md`. Everything below lives in
`claude_small_patch/` next to this file and in the data folders listed in Section 2.
The previous study of the single-column patch is `claude_column_patch_study.md` (2026-09-29, revised 2026-10-13);
this file does not repeat it.

Status markers: **[done]** finished and checked, **[running]** still computing when this line was written,
**[open]** not done.

---

## 0. Answers in short

**1. Is the single-column patch still robust in the paper's own tests (delta = C_1/dt)?** Yes, in the sense of
the paper: it converges in every run in which the star patch converges (the capped runs are the same, except
N_z = 40 at dt = 1 s in the linear dz sweep, where the column converges in 96 iterations and the star hits the cap),
and its counts stay bounded. In the linear code it reproduces the paper's figs. 4-7 almost exactly (median ratio to the
star 1.00 in all four sweeps; identical counts at dt = 50 s for every mesh, within one iteration for every aspect
ratio), with one weak spot,
thin domains L_z <= 1 km at dt = 500-1000 s (up to 2.75 x). In the nonlinear code (figs. 8-12) it equals the
star for K >= 1 (dt <= 100 s with the paper's C_1 = 1e-5) and costs two to three times the iterations (median
2.1-2.9, up to 5.9) below the threshold (dt >= 500 s). Data and the paper-style plots are in the new folders of
Section 2.2; the comparison tables and figures in `claude_small_patch/results/`.

**2. The two ways to build the smaller patch.** `LineSmoothPC` with `codims 0,1` is the single-column patch
with its two facet lines split into separate blocks (non-overlapping block Jacobi: 284 + 40 unknowns instead of
364); it has exactly the same divergence-free content, hence the same kernel condition and constants, and about
two thirds of the factorisation cost. It needs a 20-line subclass because `ASMLinesmoothPC` refuses mixed spaces.
The supervisor's "star with smaller column dofs + column Jacobi additively" is, read literally, the same operator
(a `PCCOMPOSITE` of the facet-line and column-interior blocks); read as "2D vertex stars + column Jacobi" it is
different but does not cure the thin-domain constant (10 % at N_z = 20) and costs more than the star. Recipes,
option dictionaries and the reasons are in Section 3.

**3. The dt trend of fig. 14(a).** Fig. 14 was run with delta = delta_thr(dt, eps, N) (K = 1), i.e. the shift
WAS changed with dt and L_z as the analysis prescribes; the paper's figs. 7 and 12 keep C_1 = delta dt fixed,
which is K ~ 1/dt^2 along the curve, so their decreasing counts are a decreasing K. At fixed K the counts of
both codes rise to the cross-over 1/(eps N) and fall beyond it, as the damping factor rho_01 of the analysis
does; they stay bounded, so "robust in dt given a proper shift" survives. On top of that, the nonlinear code's
counts are roughly doubled at every dt and every K tested, including dt = 1 s and K = 0.1 where the damping is
negligible, because its auxiliary-operator preconditioner applies T^{-1} A_delta^{-1} instead of A_delta^{-1}: the substitution p -> p - div(u)/delta is also
applied in the shifted pressure equation, the divergence block disappears, and the pressure back-substitution is
missing. Established with diagnostic flags (not the patch, not the inner multigrid, not the nonlinearity, not the
smoother, not the preconditioning side), checked algebraically to 1e-16, and fixed in a diagnostic copy: counts
at dt = 5000 s, K = 1 drop from 15 + 14 to 8 + 7 (star) and from 15 + 16 to 8 + 7 (column); over the whole
fig. 14(a) sweep every count roughly halves and the column/star median falls from 1.67 to 1.42. The paper's
Appendix statement that the supplied operator's off-diagonal blocks equal those of A_delta does not hold for the
code. **The solver script itself was not changed**; adopting the fix (one projection per preconditioner
application; the variant the Appendix describes, substituting only in the momentum equation, restores A_delta^{-1}
only with an exact inner solve, Section 4.4) and re-running the paper's nonlinear figures is Molin's decision.
Section 4.

**4. This file** is the record; everything created is listed in Section 1.

---

## 1. What was changed and what was created

### 1.1 Edited files (default-preserving, except that `claude_paper_run_one.sh` now deletes VTK output by default; `git diff` shows the three tracked ones, `run_claude_paper.sh` and `claude_paper_run_one.sh` are untracked (`??` in `git status`); nothing is committed)

| File | Change | Why |
|---|---|---|
| `LB_time_slice.py` | new flag `--patch {star,vanka}` (default `star`): the `mg_levels` patch block of the Schur-complement solver is `firedrake.ASMStarPC, pc_star_construct_dim 0` (as before) or `firedrake.ASMVankaPC, pc_vanka_construct_dim 1`; the `--reordering` keys use the matching prefix (`pc_star_` / `pc_vanka_`). `solver_name` prints `MG ASMVanka` for the latter. | the single-column patch existed only in the (uncommitted) nonlinear script; the paper's linear figures (figs. 4-7) need it in the linear code |
| `Nonlinear_Boussinesq_Irk_SC.py` | (the `git diff` also contains the uncommitted `--patch {star,vanka}` flag, `patch_params` block and `MG ASMVanka` solver name of the 2026-09-29 study) plus diagnostic flags, all off by default: `--linear_forms` (use `utils.LB_*_Irk`, the linear Boussinesq forms, in the equation and in the auxiliary-operator preconditioner), `--amplitude` (initial buoyancy amplitude, default 3e-4), `--mg_smoother {gmres,chebyshev,richardson}` (default gmres), `--fs1_lu` (exact LU for the velocity-buoyancy block instead of the V-cycle), `--outer_gmres_left` (outer `gmres` with left preconditioning instead of `fgmres` right). | item 3: to separate the causes of the dt trend. With no flag the solver dictionary is unchanged (checked: same keys and values; the dx sweep at dt = 50, Nx = 20 reproduces the earlier totals) |
| `utils.py` | `LB_velocity_Irk(..., use_rotation=True)`: the Coriolis term can be switched off; default unchanged | needed by `--linear_forms` with `--no_rotation` |
| `run_claude_paper.sh` | new options `-g GENERATOR` (default `Bash/bash_claude_paper.py`), `-D DATA_PREFIX` (default `data_claude_paper`), `-S SH_PREFIX` (default `test_claude_paper`); the plot command is now found by `^python .*plot_[A-Za-z_]*\.py` instead of `plot_claude_paper` only | run the paper's own sweeps (another generator, the paper's folder naming) through the same parallel driver |
| `claude_paper_run_one.sh` | deletes the VTK output of a run (`Nonlinear_slice_imp_ASM*`, `lb_slice_imp_ASM*`, ~5 MB per run) from its job directory as soon as the run has finished, unless `KEEP_VTK=1` | the disk was at 98 % (4.7 GB free) when the session started; the 2026-09-29 sweeps left 2.2 GB of it in `jobs/` |

### 1.2 New files

| Path | What |
|---|---|
| `claude_small_patch_robustness.md` | this file |
| `Bash/bash_paper_patch.py` | command generator for the paper's sweeps (delta = C_1/dt, the parameter lists of `Bash/bash.py` and `Bash/bash_nonlinear.py`) with `--patch {star,vanka}`; prints the paper's plot command (`plot_dx.py`, ...). With `--patch star` its `mpiexec` lines equal those of `test_{dt,ar,dx,dz}.sh` and `test_nonlinear_*.sh` character for character apart from the explicit solver flags listed in its docstring and a trailing `  # C1=<value>` comment on every line (checked for all eight sweeps: 0 differing `mpiexec` lines once both are stripped). |
| `claude_small_patch/run_all_sweeps.sh` | runs the four sweeps for vanka then star, nonlinear then linear (Section 2). Only the nonlinear single-column sweeps were produced through it (`run_all_sweeps.log`, 10 concurrent runs, 05:45-07:34); the linear sweeps, vanka then star, were launched separately in parallel (`run_linear.log`, 4 concurrent runs, 06:04-08:27) and the nonlinear star sweeps afterwards (`run_nonlinear_star.log`, 12 concurrent runs, 09:51-11:16); the pcfix sweeps are in `run_pcfix_ar.log` (3 concurrent runs) |
| `claude_small_patch/diag_item3/` | `run_diag.sh`, `run_diag2.sh`, `run_diag3.sh` and their logs `run_diag*.log`, `commands*.txt` (diagnostic runs of Section 4, one directory per run under `runs/`, 52 directories), `analyse_existing.py` (per-Newton-step tables and the C_1-to-K conversion from existing data; writes `per_step_counts_and_C1_vs_K.md`, `fig_C1_trajectories.png`, `fig_per_step.png`) |
| `test_paper_<mode>_<test>_<suffix>.sh` | serial forms of the new sweeps, written by the driver (same shape as the old `test_*.sh`) |
| `claude_small_patch/compare_paper_patches.py` | item 1 analysis: reads the paper's folders and the new ones, writes `claude_small_patch/results/` (`summary_<mode>_<test>.md`, `counts_*.png`, `history_*.png`, `overview_<mode>.png`, `all_counts.csv`) |
| `claude_small_patch/diag_item3/Nonlinear_pcfix.py` (+ symlink `Nonlinear_Boussinesq_Irk_SC.py` -> it) | diagnostic copy of the solver with `--pc_fix` (class `HDivSchurPCFixed`) and `--pc_fix2` (the Appendix's variant, Section 4.4); the symlink lets `run_claude_paper.sh ... -x "--script_dir ../Boussinesq_implicit/claude_small_patch/diag_item3 --extra=--pc_fix"` run it |
| `claude_small_patch/diag_item3/check_aux_factorisation.py`, `.txt` | numerical check of P_aux = A_delta T (Section 4.4) |
| `claude_small_patch/diag_item3/plot_fig14a_pcfix.py`, `fig14a_pcfix.png`, `fig14a_pcfix.md` | fig. 14(a) with and without the back-substitution |
| `claude_small_patch/item2_scripts/` | the item-2 agents' throwaway scripts (dof/patch enumeration, LU fill, 1D splitting constants) |
| `claude_small_patch/paper_style.py` | the one matplotlib style used by every figure below (fonts, viridis ramp, solid/filled = star, dashed/open = single column, figure sizes) |
| `claude_small_patch/restyle_paper_figures.py` | the ten iteration figures of the paper regenerated in that style, same file names (Section 7) |
| `claude_small_patch/restyle_fig14.py` | fig. 14 (`column_patch_sweeps.png`) in that style, same data |
| `../Patch-Preconditioner-for-Boussinesq/Image_restyled/` | the restyled figures (Section 7); the paper's `Image_new/`, `Image_nonlinear/`, `Image_patch/` are untouched |
| `../data_claude_paper_nonlinear_ar_patch{star,vanka}_pcfix/` | the K = 1 aspect-ratio sweep with the corrected preconditioner (Section 4.5) |
| `Bash/bash_claude_paper.py` | one new option `--extra` (text appended to every solver command), used for the pcfix sweeps |

### 1.3 Not touched

`Bash/bash.py`, `Bash/bash_nonlinear.py`, `plot_*.py`, `claude_paper_shift.py`, `plot_claude_paper.py`, the paper, all earlier data folders.

---

## 2. Item 1: the paper's robustness tests with the single-column patch

### 2.1 What is run

The sweeps are exactly those behind the paper's figures (not the K sweeps of `claude_paper_tests_README.md`):
the shift is delta = C_1/dt with the paper's C_1, and the lists are those of Molin's generators.

| code | test | swept | fixed | runs |
|---|---|---|---|---|
| linear `LB_time_slice.py` (`bash.py`) | dt | C_1 in {1, 0.1, 0.01, 1e-3, 1e-4, 1e-5} | nominal 20 x 20, 300 km x 4 km | 54 |
| | ar | L_z in {12, 8, 4, 2, 1, 0.5, 0.25} km | C_1 = 1e-3 | 63 |
| | dx | N_x in {10, 20, 40, 60, 80} | C_1 = 1e-3 | 45 |
| | dz | N_z in {10, 20, 40, 60, 80} | C_1 = 1e-3 | 45 |
| nonlinear `Nonlinear_Boussinesq_Irk_SC.py` (`bash_nonlinear.py`) | dt | C_1 in {1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6} | nominal | 70 |
| | ar | L_z in {12, 8, 4, 2, 1} km | C_1 = 1e-5 | 50 |
| | dx | N_x in {10, 20, 40, 60, 80} | C_1 = 1e-5 | 50 |
| | dz | N_z in {10, 20, 40, 60, 80} | C_1 = 1e-5 | 50 |

dt in {1, 5, 10, 50, 100, 500, 1000, 5000, 1e4} s (linear) and the same plus 5e4 s (nonlinear); `maxit 100`; two
time steps per run, the second one monitored; 2 MPI ranks per run; 10 runs at a time for the nonlinear single-column
sweeps, 4 for the linear sweeps and 12 for the nonlinear star re-run (the `concurrency` line of each log).

Solver flags that reproduce the paper's runs with the scripts at HEAD (see `claude_paper_tests_README.md` Section 6
for why each is needed): linear `--richardson --reordering --rtol 1e-8` (rotation off by default there);
nonlinear `--no_rotation -firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2 -firedrake_1_snes_atol 1e-9`
(f = 0, Eisenstat-Walker on as in `data_nonlinear_*_deg2_ew`, the paper's data). Plus `--patch vanka` for the
single-column patch; the star runs are repeated with the same flags as a like-for-like baseline (the paper's own
folders are the second baseline; the 2026-09-21 regression runs showed that the code at HEAD reproduces them).

### 2.2 Folders (next to the paper's `data_nonlinear_<test>_deg2_ew` and `data_linear_<test>_deg2`)

| folder | content |
|---|---|
| `../data_nonlinear_<test>_deg2_ew_vanka` | nonlinear, single-column patch **[done]** 220/220 runs |
| `../data_linear_<test>_deg2_vanka` | linear, single-column patch **[done]** 206/207 (the C_1 = 1, dt = 1 s run breaks down, see 2.3) |
| `../data_nonlinear_<test>_deg2_ew_star` | nonlinear, star patch, same flags **[done]** 219/220 (the dt = 1 s, C_1 = 1e-3 run stops with `DIVERGED_MAX_IT` after the 10 Newton steps allowed by `snes_max_it`; the paper's folder has no `snes_error_`/`snes_ksp_cum_`/`residual_` file for it either, only a stale `error_dt1.0_shift0.001.out` of 2026-05-02 with 149 entries, above the cap) |
| `../data_linear_<test>_deg2_star` | linear, star patch, same flags **[done]** 206/207 (same breakdown) |

Each folder holds the solver's `.out` files (`error_`, and for the nonlinear code `residual_`, `snes_error_`,
`snes_ksp_cum_`), the paper's plots made by `plot_dx.py` etc. (`error_dx_Robust_dt50.0.png`,
`snes_error_dx_dt50.0.png`, `error_shift.png`, ...), `commands.txt`, `run_all.log`, and `jobs/job_<hash>/` with the
command, log and status of every run (no VTK output, see Section 1.1).

### 2.3 Results

Tables with every count: `claude_small_patch/results/summary_<mode>_<test>.md` (cells "paper | star re-run |
single column"); figures `counts_<mode>_<test>.png` (counts against dt), `history_<mode>_<test>_dt<dt>.png`
(the paper's error-history figures for both patches at the paper's dt), `overview_<mode>.png` (four panels in the
style of the paper's fig. 14, for the paper's parametrisation); the paper-style plots of each folder
(`error_dx_Robust_dt50.0.png`, `snes_error_AR_dt100.0.png`, `error_shift.png`, ...) are in the data folders,
made by Molin's `plot_*.py`.

**Linear code** (the paper's figs. 4-7; C_1 = 1e-3 for ar/dx/dz; KSP iterations to rtol 1e-8):

| sweep | single column / star (runs where both converge) | not converged (star / column) | where they differ |
|---|---|---|---|
| dx | median 1.00, range 0.77-1.03 (42 runs) | 3 / 3 (N_x >= 40 at dt = 1 s, capped at 100 in both) | nowhere for dt >= 50 s: 8, 6, 5, 5, 4, 4 for every N_x with both patches |
| dz | median 1.00, range 0.88-1.33 (42) | 3 / 2 | N_z >= 40 at dt = 1000 s: 5 against 4; N_z >= 60 at dt = 10^4 s: 4 against 3 (the 1.33); equal at dt = 5000 s; one iteration either way at dt <= 10 s, and N_z = 40 at dt = 1 s converges (96) where the star is capped |
| ar | median 1.00, range 0.83-2.75 (61) | 2 / 2 (L_z >= 8 km at dt = 1 s) | **thin domains at dt = 500-1000 s**: L_z = 0.25 km 11 / 11 against 4 / 4, 0.5 km 10 / 9 against 4 / 4, 1 km 8 / 7 against 4 / 4 (star); at dt = 5000 s 6 against 3 |
| dt | median 1.00, range 0.93-1.75 (46) | 8 / 8 (the same runs: C_1 >= 0.01 at small dt, capped; C_1 = 1 at dt = 1 s breaks down, `DIVERGED_BREAKDOWN`, in both re-runs where the paper's run reached the cap) | C_1 = 1e-5 at dt = 50-500 s: 6, 7, 5 against 4, 4, 3 |

So in the linear code the single-column patch reproduces the paper's four figures: the dx and dz figures at
dt = 50 s are indistinguishable (8 iterations for every mesh with both patches, overlapping histories), the
aspect-ratio figure at dt = 50 s agrees to within one iteration (19, 13, 8, 5, 5, 5, 5 against 19, 14, 8, 6, 5, 4,
4), and the shift figure is the same except for the two smallest C_1 (C_1 = 1e-4: one iteration more at every
dt >= 100 s, 5 against 4 and 4 against 3; C_1 = 1e-5: 6, 7, 5 against 4, 4, 3 at dt = 50-500 s, 4 against 3
beyond). The place it costs most is the thin, weakly stratified corner (L_z <= 1 km, dt = 500-1000 s), where
rho_1 = m L_z^2/(pi^2 dx^2) is smallest, as the constant of the previous study predicts; the factor there is up
to 2.75. Elsewhere the excess is at most 1.75 (C_1 = 1e-5 at dt = 100 s: 7 against 4; L_z = 2 km at dt = 500 s:
6 against 4). The star re-run reproduces the paper's data in every run except the near-capped dt = 1 s runs (73
against 80 at the nominal mesh in all four sweeps; 62 against 60 at N_x = 10, 52 against 50 at N_z = 10, 40
against 41 at L_z = 2 km; breakdown instead of the cap at C_1 = 1) and two one-iteration differences at dt = 5 s
(N_z = 60: 27 against 28; C_1 = 0.01: 63 against 64).

**Nonlinear code** (the paper's figs. 8-12; C_1 = 1e-5 for ar/dx/dz; f = 0, Eisenstat-Walker on; total KSP
iterations over the Newton steps; ratios are single column / star re-run, both run with the same flags):

| sweep | single column / star (runs where both converge) | not converged (star / column) | pattern (column against star re-run) |
|---|---|---|---|
| dx | K >= 1 (dt <= 100 s): 0.96-2.14, median 1.00 (25 runs); K < 1 (dt >= 500 s): 1.50-3.90, median 2.30 (25) | 0 / 0 | dt = 50 s (the paper's fig. 8): 24, 26, 25, 25, 24 against 25, 22, 22, 22, 22; dt >= 500 s: 13-39 against 8-11 |
| dz | K >= 1: 1.00-1.79, median 1.00 (25); K < 1: 1.80-5.88, median 2.89 (25) | 0 / 0 | dt = 50 s (fig. 9): 24, 26, 27, 27, 27 against 24, 22, 22, 22, 22; dt >= 500 s: 18-47 against 8-11, worst for N_z = 60-80 |
| ar | K >= 1: 1.00-1.80, median 1.00 (23); K < 1: 1.40-4.11, median 2.44 (25) | 2 / 2 (L_z >= 8 km at dt = 1 s, capped in both) | dt = 100 s (fig. 10 left): 21, 21, 19, 22, 18 against 19, 20, 14, 14, 10; dt = 5000 s (fig. 10 right): 16, 19, 28, 24, 23 against 10, 10, 10, 10, 10 |
| dt | K >= 1: 0.98-2.00, median 1.00 (28); K < 1: 1.07-4.14, median 2.10 (31) | 10 capped + 1 Newton divergence / 11 capped (the same runs) | dt = 500 s (fig. 11): 29, 28, 23, 24, 28, 25, 29 against 28, 27, 19, 12, 10, 9, 7 |

So in the nonlinear code the single-column patch converges in every run in which the star does, with the same
set of capped runs, and its totals stay bounded (at most 47 against the star's 8-12 in the ar/dx/dz sweeps at
dt >= 500 s and 7-16 over all K < 1 runs; in the dt sweep the star itself needs up to 28 at dt = 500 s, where
C_1 >= 5e-4 keeps K >= 1); but it
costs two to three times the iterations (median 2.1-2.9, up to 5.9) wherever the paper's shift is below the
threshold (K = C_1/(dt delta_thr) < 1, i.e. dt >= 500 s for C_1 = 1e-5: the grad-div dominated regime), and
nothing above it (median 1.00 for K >= 1, i.e. dt <= 100 s, where the two patches give identical counts in 58 of
the 101 runs).

The star re-run agrees with the paper's nonlinear data to a median ratio of 0.93-0.99 (range 0.54-1.22), lower
mainly because `-firedrake_1_snes_atol 1e-9` stops Newton after two steps where the paper's runs took a third
step of 2-3 KSP iterations on a residual already below 1e-9 (e.g. 12/3 against 9/2 at dt = 1000 s); the
Eisenstat-Walker first-step tolerance occasionally changes the Newton path by more (dt = 10 s, N_x = 40: 14/2
against the paper's 26/3). The single-column figures are the ones to compare with the star re-run, which was made
with identical flags. This is the same picture as the K = 1 sweeps of the previous study,
seen along the paper's fixed-C_1 trajectories, with one difference in presentation: the paper's nonlinear counts
include the preconditioner defect of Section 4.4 for both patches.

The division of the sweeps into "equal" and "2-4x" is the division into K >~ 1 and K <~ 1 of
`per_step_counts_and_C1_vs_K.md`; the linear code shows almost no penalty at small K because its star and column
V-cycles are both good enough for the smooth Skamarock-Klemp right-hand side at rtol 1e-8, whereas the nonlinear
code's Newton steps (and the defect of Section 4.4) expose the weaker column V-cycle.

---

## 3. Item 2: two ways to build the smaller patch in Firedrake (discussion only)

Nothing was implemented in the solver and no run was made with these smoothers. Three research agents read the
Firedrake/PETSc sources (`firedrake/preconditioners/asm.py`, `base.py`, PETSc `composite.c`, `asm.c`,
`options.c`, petsc4py `libpetsc4py.pyx`) and did the mathematics; two referees each checked every claim against
the sources and the derivations (31 of 36 claims confirmed, the rest corrected below). Their throwaway scripts are
in `claude_small_patch/item2_scripts/` (dof/patch enumeration with Firedrake's own `build_*_indices`, LU fill,
the 1D splitting constants). Setting as in `claude_column_patch_study.md` Section 2: on the (u, u_y, b) block,
after eliminating u_y and b, a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w), c = dt/(2 delta) = C_M^2/K,
gamma = c/dx^2 (162 on the nominal fine mesh at K = 1); divergence-free kernel = curl of CG_2(x) x CG_2(z)
stream functions; the column patch holds only the x-bubble part of it and works for k = 2 because the coarse
level supplies the rest (constant C(rho) from 12.3 to 2, rho_1 = m L_z^2/(pi^2 dx^2)).

### 3.1 Bookkeeping: what the candidate patches contain (k = 2, N_z layers; 364 / 608 at N_z = 20)

| block | dofs | kernel functions inside |
|---|---|---|
| facet line i (base vertex): u_x normal fluxes on the vertical facet line | 2 N_z = 40 | **none** |
| cell column j, interior (base cell): u_x x-bubbles 2N_z, u_z 4N_z+2, u_y 4N_z, b 4N_z+2 | 14 N_z + 4 = 284 | curl(bubble_j(x) g(z)): (k-1)(kN_z-1) = 39, **the same as the Vanka column** |
| Vanka column (A, `ASMVankaPC construct_dim 1`) = interior + both facet lines | 18 N_z + 4 = 364 | 39 |
| star of a base vertex (committed) = facet line + two interiors | 30 N_z + 8 = 608 | (2k-1)(kN_z-1) = 117 (adds the hat-type curl(phi_i(x) g(z))) |
| 2D vertex star of the extruded mesh (`ASMExtrudedStarPC construct_dim 0`): 4 cells, 4 facets, 1 vertex, two layers tall | 52 (30 at the lids) | 9 per star; the union spans the whole kernel |

The set identities were verified with Firedrake's own index builders on a 4 x 3 slice: the codims 0,1 blocks are
disjoint and cover every dof once; Vanka(j) = interior(j) + facet(j) + facet(j+1); star(i) = facet(i) +
interior(i-1) + interior(i). Why the interior block has the same kernel as the Vanka column: a patch field with
a non-zero flux through a lateral facet has non-zero divergence (phi_i'(x) q(z), phi_i' a non-zero polynomial) in
the neighbouring column, where the patch has nothing to compensate it, so every divergence-free field of the Vanka
column already has zero lateral fluxes; the facet line alone contains no divergence-free field at all.

### 3.2 Suggestion (2): `LineSmoothPC`

`firedrake.ASMLinesmoothPC` (asm.py 255-310) loops, for each requested codimension, over the points of the BASE
DMPlex (`mesh.topology_dm` of an extruded mesh is the base interval), skips ghost points (`pyop2_ghost` label),
and takes `section.getDof(p)`/`getOffset(p)`: on an extruded mesh the section lives on base points and
`getDof(p)` counts all layers, so one patch = all dofs in the column above a base cell (codim 0) or a base
vertex (codim 1). With `pc_linesmooth_codims '0,1'` the patches are exactly the two blocks of 3.1, disjoint and
covering, so PCASM (type BASIC, asm.py 58) is **non-overlapping block Jacobi** with blocks (cell-column interior,
facet line). This is the Vanka column patch with its two facet lines split off into their own blocks.

Two practical points:

* **It refuses mixed spaces**: `if len(V) > 1: raise NotImplementedError("Not implemeneted for mixed spaces")`
  (asm.py 288-289), and `fieldsplit_1` (`pc_fieldsplit_1_fields '0,1,2'`) is a three-component MixedFunctionSpace
  on every multigrid level; line 283, a few lines earlier, even reads a section (`V.dm.getDefaultSection()`) from
  the mixed DMShell that has none. Every
  other part of `ASMPatchPC.initialize` is mixed-aware (the lgmap, the tinyasm branch, the view), so a subclass
  overriding `get_patches` is enough; `get_entity_dofs` (asm.py 582-602), the helper the working star patch uses,
  already does the per-subspace loop:

```python
import numpy
from firedrake.petsc import PETSc
from firedrake.preconditioners.asm import ASMPatchPC, get_local_ises_indices, get_entity_dofs

class MixedLinesmoothPC(ASMPatchPC):
    """ASMLinesmoothPC for mixed spaces on extruded meshes: one patch per owned base-mesh point of the
    requested codimensions (0 = base cell -> cell column interior, 1 = base vertex -> vertical facet line),
    holding every dof of every subspace in the column above that point."""
    _prefix = "pc_mixedline_"

    def get_patches(self, V):
        mesh = V.mesh().unique()
        assert mesh.cell_set._extruded
        dm = mesh.topology_dm
        codims = [int(c) for c in PETSc.Options(self.prefix).getString("codims", "0,1").split(",")]
        idx = get_local_ises_indices(V)
        ises = []
        for codim in codims:
            for p in range(*dm.getHeightStratum(codim)):
                if dm.getLabelValue("pyop2_ghost", p) != -1:
                    continue                       # owned points only: each patch is built on exactly one rank
                dofs = numpy.asarray(get_entity_dofs(V, idx, [p]), dtype=PETSc.IntType)
                dofs = dofs[dofs >= 0]
                if dofs.size:
                    ises.append(PETSc.IS().createGeneral(dofs, comm=PETSc.COMM_SELF))
        return ises
```

  and the `patch_params` block of the solver becomes (keys relative to `mg_levels`; the GMRES(6) wrapper stays)

```python
patch_params = {
    "pc_type": "python",
    "pc_python_type": __name__ + ".MixedLinesmoothPC",   # same mechanism as __name__ + '.HDivSchurPC'
    "pc_mixedline_codims": "0,1",                         # cell-column interiors (284) + facet lines (40)
    "pc_mixedline_sub_sub_pc_type": "lu",
    "pc_mixedline_sub_sub_pc_factor_mat_ordering_type": "rcm",
    "pc_mixedline_sub_sub_pc_factor_reuse_ordering": None,
    # "pc_mixedline_sub_pc_asm_local_type": "multiplicative",  # rank-local block Gauss-Seidel in the listed
    #                                                            order (cells, then facet lines); no transpose apply
}
```

  (prefix chain: `pc.getOptionsPrefix() + 'pc_mixedline_'`, the inner PCASM `..._sub_`, its block KSPs
  `..._sub_sub_`; Firedrake's defaults `sub_pc_type lu` / natural ordering are set only when absent, so `rcm` wins.
  As with every key of the `mg_levels` dictionary, these are most conveniently set in the solver dictionary; on the
  command line they work too but need the full auto-generated prefix, `-firedrake_1_...`, as `-firedrake_1_snes_atol
  1e-09` in Section 4.1.)

* **Zero-code alternative with existing classes**, composed by `PCCOMPOSITE`: the cell-column interiors are
  exactly `firedrake.ASMStarPC` with `pc_star_construct_dim 1` (the star of a base cell is the cell itself: 14N_z+4
  mixed dofs, verified), and the facet lines are `firedrake.ASMVankaPC` with `pc_vanka_construct_dim 0`,
  `pc_vanka_exclude_subspaces '0,1,2'`, `pc_vanka_include_type entity` (every subspace restricted to the seed
  vertex: 2N_z dofs, verified):

```python
patch_params = {
    "pc_type": "composite", "pc_composite_type": "additive", "pc_composite_pcs": "python,python",
    "sub_0": {"pc_python_type": "firedrake.ASMStarPC", "pc_star_construct_dim": 1,          # cell columns
              "pc_star_sub_sub_pc_type": "lu", "pc_star_sub_sub_pc_factor_mat_ordering_type": "rcm",
              "pc_star_sub_sub_pc_factor_reuse_ordering": None},
    "sub_1": {"pc_python_type": "firedrake.ASMVankaPC", "pc_vanka_construct_dim": 0,        # facet lines
              "pc_vanka_exclude_subspaces": "0,1,2", "pc_vanka_include_type": "entity",
              "pc_vanka_sub_sub_pc_type": "lu"},                                               # (diagonal block)
}
```

  Mechanics (checked in `composite.c`): `PCSetUp_Composite` hands the parent's DM and operators to every child
  (lines 165-178), so Firedrake python PCs, which read `pc.getDM()` for the function space, work as children;
  `PCSetFromOptions_Composite` (239-267) creates the children of `pc_composite_pcs` through `PCCompositeAddPCType`,
  whose helper gives them the prefixes `sub_0_`, `sub_1_` (372-375), and then calls `PCSetFromOptions` on every
  child (261-265); the additive apply is y = sum_i P_i x (125-146), `multiplicative` is
  y = P_1 x + P_2 (x - A P_1 x) in the listed order, `symmetric_multiplicative` adds the reverse pass; there are no
  weights. One trap: `PCSetFromOptions_Composite` appends new children on every call, so it must be called exactly
  once, which `PCSetUp_MG` guarantees for the level smoothers (`if (!pc->setupcalled)`), but a second
  `setFromOptions` (a reset) would duplicate them. Firedrake's regression tests
  (`tests/firedrake/regression/test_star_pc.py`) do exercise both classes on mixed spaces (`ASMStarPC` on a
  `MixedFunctionSpace`, `ASMVankaPC` with `pc_vanka_exclude_subspaces`), but not on an extruded mesh with
  `construct_dim 1`, nor with `include_type entity` applied to every subspace, so this configuration rests on
  reading the code, not on a shipped test.

Parallel: patches are built on owned base points only (the `pyop2_ghost` test), each column or facet line on
exactly one rank, no empty index sets (the star/Vanka builders emit one per ghost seed), and all dofs of an owned
point are owned, so the PCASM scatters are local. The `(VERTEX, 1)` overlap of the script is still needed for the
star patch and for assembling the interior-facet terms.

### 3.3 Suggestion (1): "a star patch with smaller column dofs, and Jacobi iteration outside on the column additively"

Three readings were examined.

**R1 (the reading under which (2) "performs this").** Keep of the base-vertex star only the dofs on the vertex's
own column, the facet line (2N_z fluxes), and treat the remaining dofs, the cell-column interiors, by block Jacobi
added additively. "Outside" is then the complement of the reduced star, the two families partition the dofs, and
the operator is sum_i R_i^T A_i^{-1} R_i over facet lines and interiors: **exactly `ASMLinesmoothPC codims 0,1`
of 3.2** (PCASM-BASIC on disjoint blocks, additive composition = union of the families; the composite above is the
literal construction). Ways (1) and (2) coincide; the composite only adds the multiplicative orderings (one block
Gauss-Seidel sweep: fluxes, then exact column solves on the updated residual) or, with a `pcksp`-Richardson child,
a weight.

**R2 ("smaller" = fewer layers).** The genuine 2D vertex star of the extruded mesh, `firedrake.ASMExtrudedStarPC`
with `pc_star_construct_dim 0` (asm.py 381-519: base vertex star x 1D vertex star = 1 vertex, 2 vertical and 2
horizontal facets, 4 cells; 52 dofs; N_x(N_z+1) patches per level, 1680 on the finest; it loops over the subspaces,
so it works on the mixed block), composed additively with the column blocks of 3.2 to supply the vertical coupling
the small patch lacks (recipe: `pc_composite_pcs python,python`, `sub_0` = `firedrake.ASMExtrudedStarPC`,
`pc_star_construct_dim 0`, optionally `pc_star_use_coloring True` (4 block-diagonal index sets if every level has an
even N_x, as here); `sub_1` = the cell-column blocks). This is genuinely different from (A)/(B): the vertex stars
hold the hat-type kernel functions the column lacks, in two-layer pieces, and on their own satisfy the kernel
condition for every k without a coarse space (Arnold-Falk-Winther) - robust in c and K.

**R3 (the most literal reading; conjecture).** The star's advantage is the hat-column functions
curl(phi_i(x) g(z)) = (phi_i g', -phi_i' g), which need only the facet line i (u_x) and the u_z (and b) dofs of the
two adjacent columns: facet line + u_z, b of two columns = 18N_z+8 = 368 dofs, plus cell-column Jacobi for the
bubbles. The entity-based ASM classes cannot separate u_x-interior from u_z dofs (same base-cell point in the
enriched RT element); it would need either a component-aware `get_patches` or the RT element split into its two
HDiv parts as separate fields (then `ASMVankaPC construct_dim 0`, `exclude_subspaces` = {u_h, u_y}, `include_type
entity`). Whether it keeps the star's constant at finite c is argued (same nodal splitting in x available), not
shown.

"Jacobi" has to mean **block** Jacobi over whole columns. Point Jacobi on the column interior fails: every single
RT basis function has non-zero divergence, so the bubble kernel modes curl(beta_j g), which no coarse correction
removes, cost c/((m + ...) dz^2) ~ 3e4 at K = 1, dt = 50 s on the nominal mesh in the dof-wise decomposition
(worse than the one-level column patch's 1/(15 gamma)), and the field-1 solve stalls. The only place a diagonal
treatment is harmless is where R1 actually uses it, on the facet lines, which carry no kernel: there the facet
block is a scalar multiple of the DG_1(z) mass matrix (diagonal with Firedrake's Gauss-Legendre DG_1 nodes) for the
rest-state linearisation, and diagonal to 1e-9 with the nonlinear advection terms at amplitude 3e-4.

### 3.4 Why (B) does what is wanted, and what it does not

1. **Kernel condition.** N_h ∩ interior_j = N_h ∩ Vanka_j and N_h ∩ facet_i = {0} (3.1), so (B) satisfies the
   Schoeberl / Lee-Wu-Xu-Zikatanov condition N_h = N_2h + sum (N_h ∩ V_i) exactly when (A) does: for k = 2 through
   coarse level + x-bubbles, not for k = 1. The two-level constant of the kernel splitting, C(rho) = 12.28, 11.71,
   8.50, 3.35, 2.19, 2.02, 2.00 at rho = 0, 1e-3, 1e-2, 0.12, 1, 10, infinity, and its growth with levels (12.3,
   15.7, 36, 74), carry over to (B) number for number (an agent's independent 1D model reproduces the study's
   table to every digit). So (B) inherits the column patch's robustness profile: bounded as K -> 0 for k = 2,
   worst (12-16) in thin, weakly stratified domains, growing with the number of levels there, nothing for k = 1.
2. **The complement.** The interiors cannot change a column's total divergence (div of the interior block =
   zero-mean functions in the column), so the column means are a 1D problem in the depth-integrated fluxes F_j
   with energy c ||D_x F||^2. (B) relaxes it with Jacobi on single facet lines (smoothing factor 1/3 with optimal
   damping), (A) with overlapping pairs {F_j, F_j+1} (0.059); x-smooth flux error is left to the coarse grid in
   both. Under the GMRES(6) wrapper this difference should be modest (within the +-20 % run-to-run noise the
   study reports for column counts at dt <= 100 s); it is the one thing a run would have to confirm. Both have
   chromatic number 3 (facet lines i and i+1 are coupled through the column between them), so lambda_max is
   bounded alike.
3. **Cost.** One banded LU of 14N_z+4 per column plus a diagonal facet block, against one of 18N_z+4 whose 4N_z
   facet dofs are factorised twice. With the sparsity pattern of the solver's Jacobian (DG upwind terms: PETSc
   reports 39 632 stored entries in the 608-dof star block and 53 468 in its LU, `claude_column_patch_study.md`
   Section 3), a scipy/RCM factorisation of the same pattern gives fills of about star 64.8k, Vanka 23.2k, (B)
   15.3k per base column; the two factorisations count differently (64.8k against PETSc's 53.5k for the same
   block), the ratios are what matters: (B) is 0.66 of Vanka and 0.24 of the star. (The saved
   `item2_scripts/lu_fill.py` assembles the rest-state form without advection, i.e. the earlier 13.7k / 5.8k
   setting; the solver-pattern numbers are the referee's re-run of it with the upwind terms, not kept.) Since `MatSolve` was 0.056 s of the column patch's 0.21 s per FGMRES
   iteration, the saving is 5-11 % per iteration, partly eaten by twice as many PCASM subdomains (2N_x).
4. **Prediction for the sweeps.** (B) within noise of (A): the same failures, the same 1.0-1.7 median ratios to
   the star, the same 3-10 in thin domains for the Skamarock-Klemp data. It does not remove the column patch's
   drawbacks; it is a cheaper way to the same operator.

### 3.5 Why R2 does less than one might hope, and what would restore the star's constant

The vertex stars alone are robust in c but not in the vertical anisotropy: on the kernel the operator is
(psi_z, psi_z) + m (psi_x, psi_x), with cell anisotropy eta = m dz^2/dx^2 = 2.8e-3 m on the nominal fine mesh,
and with coarsening in x only no coarse space ever acts in z. A hat-column kernel function curl(phi_i(x) g(z))
smooth in z must be split layer by layer into two-layer star pieces with no coarse space in z, at the price
C_s(N_z, eta) ~ N_z^2/(pi^2 (1 + pi^2 rho_1)) (exact values of the 1D model in `item2_scripts/asm_constants.py`:
81 = 2N_z^2/pi^2 at N_z = 20 for eta -> 0, twice the asymptotic formula; 26 at the nominal dt = 50 s (eta = 3.0e-3);
~410 at N_z = 80 at the same rho_1 = 0.12; O(1) once eta >~ 0.1, i.e. dt >~ 1000 s on the nominal geometry), times
an x-factor 1.5. Adding the column blocks of (B) restores exact z-line relaxation for the complement and the bubble
route for the kernel, so each hat-type mode takes the cheaper of the two routes, 1/C_(C) >= 1/C_hier + 1/(1.5 C_s):
in the thin, weakly stratified limit, where the column patch is worst (12.3), the exact model constant of (C) is
8.4 / 11.0 / 11.9 at N_z = 10 / 20 / 40 (31 % / 10 % / 3 % better), i.e. the help vanishes with N_z, in the
direction where help is needed; it is 20-100 % for eta between 0.01 and 0.1 and substantial for vertically rough
data (random right-hand sides). Cost: the dense 52 x 52 star solves alone (N_x(N_z+1) of them, 21 x 52^2 = 57k
multiply-adds per base column) are already comparable to the star patch's MatSolve work (53-65k); with the column
blocks (15.3k) (C) is about 72k per base column per sweep, 1.1-1.4 x the star's and 3.1 x (A)'s, with 23 x the
subdomains. R2 is therefore dominated by the star on robustness and by (A)/(B) on cost; the ratios of 3-10 in
thin domains would essentially persist.

What the thin-domain constant needs is hat-type kernel functions on FULL vertical lines: either the base-vertex
star itself (608), or a potential-space line correction curl(phi_i(x) x CG_2^0(z)), 2N_z-1 = 39 unknowns per
vertical line, added to (B) (nodal + coarse splitting in x: model constant <= 1.9 uniformly in rho, against 12.3
for the hierarchical splitting of the column patch). That is item 5 of the previous study's conclusions
(Hiptmair-type correction) and was not asked for here.

### 3.6 Summary for the supervisor's question

| | (A) Vanka column | (B) LineSmoothPC codims 0,1 = way (1) R1 | (C) 2D stars + column blocks = way (1) R2 |
|---|---|---|---|
| kernel content | bubbles, 39/column | identical to (A) | (A) + hat pieces two layers tall |
| constant, rho -> 0 | 12.3 (2 levels), 15.7 (3) | identical | 11.0 at N_z = 20 (8.4 at 10, 11.9 at 40) |
| complement | overlapping flux pairs | Jacobi on facet lines, z-line exact | as (B) + per-layer stars |
| patch dofs per column (sum) | 364 | 324 | 1372 |
| LU fill per column (solver pattern) | 23.2k | 15.3k | ~72k (stars 57k + column blocks 15.3k), above the star |
| Firedrake | exists | 20-line subclass, or ASMStarPC dim 1 + ASMVankaPC entity hack in PCCOMPOSITE | ASMExtrudedStarPC + column blocks in PCCOMPOSITE |
| predicted counts | reference | within noise of (A), same failures | thin SK domains ~ (A); rough data better; never the star |

Both suggestions can be built with the pieces above in an afternoon; (B) is the one worth running (it is the
same operator as the Vanka column at two thirds of its factorisation cost), and its measurement would also settle
the one open point of 3.4 (the Jacobi-versus-pairs flux relaxation).

Corrections the referees made to the agents' drafts and that are incorporated above: the LU fills were first
quoted for the rest-state linear form (star 13.7k / Vanka 5.8k) and are now for the solver's denser Jacobian
pattern; the vertex-star anisotropy constants (26 and ~410 from the 1D model at the nominal eta; the asymptotic
formula's 18 and ~290 understate them, the drafts' 20-40 and 530 were wrong); (B)'s block graph is not bipartite
(chromatic number 3); the (C) cost (1.1-1.4 x star, not 0.8); the (C) constants with the x-factor 1.5
(10 % not 13 % at N_z = 20); dropping u_y from the star patch is not free (the smoother would never act on u_y);
`ASMLinesmoothPC`'s tinyasm backend stores an explicit dense inverse per patch, not an LU.

---

## 4. Item 3: the dt trend of fig. 14(a)

### 4.1 What fig. 14(a) was run with

Fig. 14 = `Image_patch/column_patch_sweeps.png`, made by `claude_column_patch/paper_figure.py` from
`../data_claude_paper_nonlinear_ar_patch{star,vanka}` (panel (a)). Every run of panel (a) is

```
mpiexec -n 2 python ../Boussinesq_implicit/Nonlinear_Boussinesq_Irk_SC.py --nx 20 --nz 20 --length 300000.0 \
    --height <Lz> --dt <dt> --tmax <2 dt> --shift <delta_thr(dt, Lz)> --maxit 100 --ar_test --no_rotation \
    [--patch vanka] -firedrake_1_snes_atol 1e-09
```

with L_z in {12, 8, 4, 2, 1, 0.5, 0.25} km, dt in {1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5} s, f = 0,
**no Eisenstat-Walker**, FGMRES rtol 1e-6 / atol 1e-9, Newton atol 1e-9, degree 2, three MG levels with six GMRES
iterations per level, and **the shift delta = K delta_thr(dt, eps, N) with K = 1** (`claude_paper_shift.py`:
delta_thr = dt/(2 C_M^2), C_M^2 = max(C_h^2, m C_v^2), C_h = L_x/(2 pi), C_v = L_z/pi, m = 1 + dt^2 N^2/4).
The quantity plotted is the total number of FGMRES iterations over the Newton steps of the monitored time step.

So yes: the shift was changed with dt (and with L_z) in fig. 14, following the analysis, whereas the paper's
figs. 7 and 12 keep C_1 = delta dt fixed.

### 4.2 Fixed C_1 against fixed K

K = delta/delta_thr = C_1/(dt delta_thr(dt)). On the nominal geometry delta_thr = 2 pi^2 dt/L_x^2 for dt below the
cross-over 1/(eps N) = 7500 s, so K = C_1 L_x^2/(2 pi^2 dt^2) **falls like 1/dt^2** along a fixed-C_1 curve
(table in `claude_small_patch/diag_item3/per_step_counts_and_C1_vs_K.md`): for the paper's nonlinear C_1 = 1e-5,
K = 4.6e4 at dt = 1 s, 18 at 50 s, 0.18 at 500 s, 1.8e-3 at 5000 s, 8e-4 beyond. For the linear C_1 = 1e-3,
K = 4.6e6 at dt = 1 s, 1.8e3 at 50 s, 18 at 500 s, 0.18 at 5000 s.

Consequence: the decrease of the counts with dt in figs. 7 and 12 is the decrease of K along the curve (the K
sweeps show the counts fall with K at every dt), not dt-robustness at a fixed distance from the threshold. At fixed
K = 1 the counts rise towards the cross-over and fall beyond it, in both codes, which is the dt dependence of the
damping factor of the gravest mode present, rho_01 = K/(K + 2/min(1, 4 eps^2 m)) (3.6e-4 at dt = 1 s, 0.18 at
5000 s, 1/3 from the cross-over on). The two statements are consistent: "robust in dt given a proper shift" holds
in the sense that at K <= 1 the count is bounded for all dt (linear code 3-5; nonlinear 2-15 per Newton step);
the shape of the curve depends on how the shift is tied to dt.

`fig_C1_trajectories.png`: the iteration maps of the two K sweeps with the paper's fixed-C_1 lines overlaid.

### 4.3 Where the nonlinear code's extra iterations come from

Per-Newton-step counts (star patch, nominal geometry, K = 1, no EW; `per_step_counts_and_C1_vs_K.md`):

| dt | 1 | 10 | 50 | 100 | 500 | 1000 | 5000 | 1e4 | 5e4 | 1e5 |
|---|---|---|---|---|---|---|---|---|---|---|
| Newton step 1 | 4 | 4 | 4 | 5 | 6 | 9 | 15 | 13 | 9 | 8 |
| Newton step 2 | - | 2 | 9 | 7 | 11 | 11 | 14 | 15 | 10 | 10 |
| linear code `LB_time_slice.py`, K = 1, rtol 1e-8 | 4 | 3 | 3 | 4 | 4 | 4 | 5 | 5 | 3 | 3 |

In the aspect-ratio sweep of fig. 14(a) the first-step count of the star patch behaves the same way for every
L_z, with the maximum moving to the cross-over 1/(eps N) of that L_z (e.g. L_z = 250 m: 3, 3, 3, 3, 4, 4, 8, 9,
14, 16 for dt = 1 ... 1e5 s, cross-over 1.2e5 s), i.e. the shape is the rho_01(dt, eps) shape of the analysis.

Diagnostic runs at dt = 5000 s, K = 1, star patch, nominal geometry (`claude_small_patch/diag_item3/runs/`;
numbers are KSP iterations of Newton step 1 + step 2 [+ step 3]):

| variant | counts | what it rules out |
|---|---|---|
| committed solver (nonlinear forms, V-cycle inner solve, FGMRES right, rtol 1e-6) | 15 + 14 | (the sweep value, 29) |
| `--linear_forms` (linear equations, everything else the same) | 13 + 5 | **not the nonlinearity**: with linear forms the first step is 4, 4, 6, 13, 13, 7 at dt = 1, 50, 500, 5000, 1e4, 1e5, the same rise as with the nonlinear forms |
| `--linear_forms --fs1_lu` (exact LU inner solve) | 13 + 2 | **not the inner multigrid** (the same 13-iteration staircase; the residuals agree to 1e-6 relative over the first ten iterations, the converged residual is 2.6e-7 against 2.7e-8) |
| `--fs1_lu` (nonlinear forms) | 15 + 13 | idem |
| `--linear_forms --mg_smoother chebyshev` (the linear code's smoother) | 14 + 3 | not the smoother |
| `--amplitude 3e-5` / `3e-3` (nonlinear forms) | 13 + 10 / 15 + 14 + 11 | first step independent of the data amplitude |
| `--linear_forms`, K = 0.5 / K = 0.1 | 11 + 5 / 8 + 4 | weak dependence on K (K_eff = 2K, derivation at the end of Section 4.4) |
| `--linear_forms --rtol 1e-8` | 16 | the last three iterations gain three decades (2.6e-7 to 2.7e-10; fast final convergence) |
| `--linear_forms --outer_gmres_left --atol 1e-13` (the linear code's outer solver) | 11 + 10 | not the preconditioning side; GMRES-left shows the staircase exactly (preconditioned-residual ratios 0.03, 1.00, 0.10, 1.00, 0.10, 1.00, 0.11, 1.00, ...; the true residual rises by 1.79 at the first iteration and then follows 1.00, 0.10, 1.00, 0.10, ...) |
| `--linear_forms --fs1_lu --outer_gmres_left --atol 1e-13` | 11 + 8 | idem with the exact inner solve |
| `--outer_gmres_left --atol 1e-13` (nonlinear forms) | 11 + 12 + 7 | idem |

The residual history of the first Newton step is the same staircase in all of them: with FGMRES the ratios are
0.92, 0.84, 0.22, 0.95, 0.10, 0.95, 0.11, 0.95, 0.07, 0.91, 0.03 (0.06 with the nonlinear forms), ..., with
GMRES-left (preconditioned residual) they are 0.03, 1.00, 0.10, 1.00, 0.10, 1.00, ...: **every second Krylov iteration achieves nothing**. The linear code at the same dt and
K converges in 5 iterations with a plain geometric history (1, 0.096, 3.8e-3, 7.8e-5, 9.3e-7, 1.5e-13). The
productive iterations of the nonlinear script (about six at dt = 5000 s) match the linear code's count.

### 4.4 The cause: the auxiliary-operator preconditioner is T^{-1} A_delta^{-1}, not A_delta^{-1}

`HDivSchurPC.getNewForm` builds the auxiliary form by the substitution p -> p - div(u)/delta **everywhere the
pressure appears, including the added shifted pressure equation** (`p = p - 1/delta * div(u)` is executed before
`F += delta * p * phi * dx`). In the pressure row the divergence then cancels, delta p phi - div(u) phi + div(u) phi
= delta p phi. Notation (stage variables of the implicit midpoint rule, a = dt/2; (pressure, velocity-buoyancy)
block ordering as in the paper's Appendix): D_div = matrix of phi div(u) dx, M_p = pressure mass matrix, A_vv =
velocity-buoyancy block (mass + a times the linearised transport and buoyancy coupling). The shifted operator the
analysis is about is

    A_delta = [[ a delta M_p ,  a D_div ],
               [ -a D_div^T  ,  A_vv    ]],

and the operator PETSc is given is

    P_aux = [[ a delta M_p ,  0 ],
             [ -a D_div^T  ,  S ]],      S = A_vv + (a/delta) D_div^T M_p^{-1} D_div = the exact Schur complement of A_delta,

lower triangular: the (pressure, velocity) block a D_div has disappeared. The two are related by a change of
variables, P_aux = A_delta T with T = [[I, -(1/delta) M_p^{-1} D_div], [0, I]] (k_p = k_p,var - M_p^{-1} D_div k_u /
delta; checked by multiplying out). The paper's Appendix statement that the off-diagonal blocks of the supplied
operator agree with those of A_delta is therefore not what the code does. The fieldsplit "full" factorisation of
P_aux does give the exact Schur complement S in its (1,1) block (Prop. aux-form is right about that), but the
operator it applies is

    P_aux^{-1} = T^{-1} A_delta^{-1},   T^{-1} = [[I, (1/delta) M_p^{-1} D_div], [0, I]],

i.e. the velocity of the shifted solve, but a pressure k_p,var = k_p + M_p^{-1} D_div k_u / delta = (a delta M_p)^{-1}
r_p, the pressure residual over the shift, independent of the velocity solve (the missing back-substitution of the
pressure). The effect on the preconditioned operator: with the exact Schur solve, A_delta^{-1} J = [[Lambda, 0],
[-a S^{-1} D_div^T, I]] with Lambda carrying the eigenvalues 1 - rho_jl of the analysis, whereas T^{-1} A_delta^{-1} J
= [[0, (1/delta) M_p^{-1} D_div], [-a S^{-1} D_div^T, I]] has the eigenvalue 1 on the divergence-free velocities only
and, for every rho_jl, the two roots of lambda^2 - lambda + (1 - rho_jl) = 0: a conjugate pair with real part 1/2
(e^{+-i pi/3} for rho -> 0). Every eigenvalue cluster of the analysed operator becomes a pair, GMRES needs two
iterations per cluster whatever the damping - the staircase - and the factor of about two is independent of K and
dt, as the runs below show (dt = 1 s: 4 -> 2; dt = 5000 s, K = 0.1 / 1 / 10: 16 -> 8, 29 -> 15, 58 -> 30). The
linear code (`LB_time_slice.py`) applies the fieldsplit factorisation to the natural shifted Jacobian Jp with both
off-diagonal blocks present and the grad-div form supplied separately as the Schur-complement preconditioner, so it
applies A_delta^{-1} and converges as the analysis says.

The Appendix describes F_hat with "the shifted divergence equation retained as in A_delta", i.e. the substitution
applied in the momentum equation only; that operator, P' = [[a delta M_p, a D_div], [-a D_div^T, S]], has the
off-diagonal blocks of A_delta, and the fieldsplit (full, a11, preonly on fieldsplit_1) applied to it is A_delta^{-1}
with the (1,1)-block solver in place of S^{-1}. The code deviates from this description, and the deviation is not
harmless: tried as `--pc_fix2` in the diagnostic copy, P' gives 4 + 8 at dt = 50 s and 12 + 16 at dt = 5000 s (K = 1,
star) with the V-cycle, i.e. the staircase persists (ratios 0.75, 0.86, 0.08, 0.99, 0.06, 0.99, 0.09, ...), but
8 + 7 with the exact LU inner solve (`--pc_fix2 --fs1_lu`), the same as the projection fix. So the Appendix's variant
restores A_delta^{-1} only when the Schur block is solved exactly; with the multigrid V-cycle on the Schur block the
two formulations differ in a way I did not resolve (presumably through the operator PETSc hands the level solvers
when the Schur complement of P' is not its (1,1) block). The projection fix keeps the triangular P_aux, for which the
Schur complement is exactly its (1,1) block, and works with the V-cycle; it is the one to use.

**Numerical check of the algebra** (`claude_small_patch/diag_item3/check_aux_factorisation.py`, output
`check_aux_factorisation.txt`; an agent wrote and ran it, I read the output): on a 4 x 3 mesh, degree 2, with
the stage Jacobians assembled by Irksome's `getForm` for GaussLegendre(1), ||P_aux - A_delta T||/||A_delta|| =
1e-16 for (delta, dt) = (1e-6, 50) and (1e-9, 5000); the (pressure, velocity) block of P_aux is zero (1e-16
relative); P_aux - A_delta is exactly (dt/2delta) D^T M_p^{-1} D in the velocity block (1.6e-15); the pressure
component of P_aux^{-1} r equals ((dt/2) delta)^{-1} M_p^{-1} r_p to 3e-16 and does not change when r_u changes,
while that of A_delta^{-1} r changes by 97-150 %; ||P_aux^{-1} - T^{-1} A_delta^{-1}|| = 2e-15 / 3e-13 relative.
(One bookkeeping point: with Firedrake's symmetric elimination of the Dirichlet rows and columns, T must be built
with D restricted to the free velocity columns; restricted to the free unknowns the identity holds to 1e-16.)

**Test.** `claude_small_patch/diag_item3/Nonlinear_pcfix.py` is a copy of the solver script with a flag `--pc_fix`
whose class `HDivSchurPCFixed` applies T after the auxiliary solve: k_p <- k_p - Pi(div k_u)/delta (Pi = L^2
projection onto the pressure space; the dt/2 of the stage substitution cancels), which turns the applied operator
into A_delta^{-1} with the V-cycle in place of S^{-1}. Nominal geometry, star patch, K = 1, nonlinear forms, no EW
(the setting of fig. 14), KSP iterations per Newton step:

| dt | 1 | 50 | 500 | 5000 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|
| committed preconditioner (sweep) | 4 | 4 + 9 | 6 + 11 | 15 + 14 | 13 + 15 | 8 + 10 |
| `--pc_fix` | 2 | 2 + 2 | 3 + 6 | 8 + 7 | 7 + 8 | 5 + 5 |
| linear code, K = 1, rtol 1e-8 (for comparison; by K_eff = 2K the like-for-like linear run is K = 2, not run; at K = 10 the linear code needs 4, 4, 4, 9, 8, 4) | 4 | 3 | 4 | 5 | 5 | 3 |

(the copy without `--pc_fix` reproduces the sweep's 15 + 14 at dt = 5000 s.) With the fix the residual histories are plain geometric (dt = 5000 s: ratios 0.87, 0.19, 0.10, 0.11, 0.06, 0.06,
0.30, 0.03; no idle iterations), and the other variants at dt = 5000 s, K = 1 become: linear forms 7 + 1; exact inner
LU 8 + 7 (so the V-cycle is as good as exact there); K = 0.1: 4 + 4 (was 8 + 8); K = 10: 15 + 15 (was 29 + 29);
**single-column patch 8 + 7 (was 15 + 16) at dt = 5000 s and 5 + 4 (was 14 + 10) at dt = 50 s**, against the star's
8 + 7 and 2 + 2. The two patches coincide at dt = 5000 s once the preconditioner is right; at dt = 50 s the column
patch still costs about twice the star's count, as the previous study's constant predicts. Whether the fix changes the
star/column picture of fig. 14(a) as a whole changes is checked with the full K = 1 aspect-ratio sweep for both
patches (`../data_claude_paper_nonlinear_ar_patch{star,vanka}_pcfix`, Section 4.5).

**K_eff = 2K.** With the same shift delta the linear code's preconditioning matrix (unknowns u^{n+1}, p^{n+1}) is
A_delta,lin = [[delta M_p, D_div], [-dt D_div^T, A_vv]]: its pressure row `phi div(u) + delta p phi` carries no dt,
its momentum row `-dt div(w) p` carries dt. The stage Jacobian of the nonlinear script is A_delta,nl = [[a delta M_p,
a D_div], [-a D_div^T, A_vv]], a = dt/2 (Irksome's `getForm` replaces every field by u^n + dt a k, GaussLegendre(1)
has a = 1/2), with the same A_vv because every coupling term carries dt/2 in both codes (u^{n+1/2} = (u^n + u^{n+1})/2
against u^n + a k). The pressure Schur complements are delta M_p + dt D_div A_vv^{-1} D_div^T (linear) and
a (delta M_p + a D_div A_vv^{-1} D_div^T) (nonlinear): relative to the term the shift competes with, the nonlinear
script's delta is worth 2 delta of the linear code, K_eff = 2K. Not separately measured; it is a factor of two in
K, which the K sweeps show to matter little at K <= 1.

### 4.5 Fig. 14(a) with the corrected preconditioner

The K = 1 aspect-ratio sweep of fig. 14(a) was re-run with `--pc_fix` for both patches
(`../data_claude_paper_nonlinear_ar_patch{star,vanka}_pcfix`, 70 + 70 runs, same commands otherwise; table
`claude_small_patch/diag_item3/fig14a_pcfix.md`, figure `fig14a_pcfix.png`):

| | committed preconditioner (fig. 14(a)) | with the back-substitution |
|---|---|---|
| star, totals over the sweep | 3 to 40 | 2 to 21 |
| single column | 11 to 44 | 4 to 21 |
| fixed / committed, star | | median 0.52 (0.30-0.67) |
| fixed / committed, column | | median 0.47 (0.29-0.70) |
| column / star | median 1.67, max 10.7 | median 1.42, max 7.0 |

The shape of the curves is unchanged (rise to the cross-over of each L_z, fall beyond it), every count is roughly
halved, and the star/column picture of the paper's Section 5 survives with a smaller spread: at dt = 50 s the
column patch needs 8, 8, 9, 10, 12, 13, 12 against the star's 6, 5, 4, 3, 2, 2, 2 (L_z = 12 km to 0.25 km), and
for dt >= 5000 s the two patches agree to within two iterations except in the thinnest domains.

### 4.6 Answer to the question

1. Fig. 14(a) used delta = delta_thr(dt, eps, N) (K = 1), i.e. the shift was changed with dt and with L_z as the
   analysis prescribes; the paper's figs. 7 and 12 keep C_1 = delta dt fixed, which in the analysis' terms is a
   shift that falls like 1/dt^2 relative to the threshold. The two curves cannot have the same shape.
2. At fixed K the counts of both codes rise towards the cross-over and fall beyond it; this is the predicted
   dt dependence of rho_01 at fixed K, not a loss of robustness: the counts stay bounded (linear code 3-5 on the
   nominal geometry and 3-6 over its aspect-ratio sweep, nonlinear code <= 16 per Newton step at K = 1 for every dt
   and every L_z of the sweep). "Robust in dt given a proper shift"
   remains true in that sense.
3. The nonlinear code's counts are, in addition, roughly doubled relative to the corrected preconditioner at every
   dt and K tested (dt = 5000 s, K = 0.1 / 1 / 10: 16 -> 8, 29 -> 15, 58 -> 30; dt = 1 s: 4 -> 2), because the
   auxiliary-operator preconditioner omits the pressure back-substitution (Section 4.4). This affects every nonlinear count in the paper and in the two
   studies (figs. 8-12, 14, the K sweeps of `claude_paper_tests_README.md`, and the star/column comparison of
   `claude_column_patch_study.md`), by a factor of about two that, in the runs made, does not depend on K or dt (each
   eigenvalue cluster of the analysed operator becomes a pair, Section 4.4). The fix is one projection per
   preconditioner application (the `HDivSchurPCFixed` class); whether to
   adopt it and re-run the nonlinear figures is Molin's decision - **nothing in the repository's solver was changed**.

---

## 5. Housekeeping

* Disk: `/` was at 98 % (4.7 GB free) at the start and is at 99 % (4.4 GB free) at the end. The new runs do not keep VTK output. The 2.2 GB left by the
  2026-09-29 sweeps can be freed with
  `rm -rf ../data_claude_paper_nonlinear_*_patch*/jobs/*/Nonlinear_slice_imp_ASM*` (nothing reads those files);
  I did not delete them.

---

## 6. What is solid and what is not

| Statement | Basis |
|---|---|
| Item 1 counts (Section 2.3, `results/`) | the runs themselves; the generator reproduces the old command lines character for character; the star re-run reproduces the paper's data except near the iteration cap |
| "the column patch converges wherever the star does" | every folder: identical sets of capped / broken-down runs, except linear dz N_z = 40 at dt = 1 s (star capped, column 96) |
| Item 2, Firedrake mechanics (what each class builds, mixed-space refusal, composite behaviour, prefixes) | read from the venv's sources by two agents and checked by two referees against the same sources; dof counts and set identities verified with Firedrake's own index builders on a 4 x 3 slice |
| Item 2, constants and predictions for (B) and (C) | 1D model computations (reproduce the previous study's C(rho) table to every digit) and stable-decomposition arguments; **no run of either smoother** |
| Item 3, parameters of fig. 14(a) and the C_1-to-K conversion | the command files and `claude_paper_shift.py` |
| Item 3, the cause (missing back-substitution) | diagnostic runs excluding every other candidate; algebraic identity checked to 1e-16 on the assembled stage Jacobians; the corrected copy removes the staircase and halves the counts in 140 + 13 runs |
| K_eff = 2K between the two codes | derivation from the stage scaling at the end of Section 4.4; not separately measured |
| the Appendix's variant of the fix (`--pc_fix2`) | two runs with the V-cycle (staircase persists) and one with the exact inner solve (8 + 7); the mechanism of the difference is not resolved |
| "robust in dt given a proper shift survives" | the K = 1 sweeps of both codes stay bounded; this is an interpretation of what the question meant by robust |

Not done: no run with LineSmoothPC or the composite smoothers (item 2 was discussion only); the back-substitution
fix was not applied to the solver script and no paper figure was regenerated with it. The 2.2 GB of VTK leftovers of the
2026-09-29 sweeps were not deleted (Section 5).

---

## 7. Figures for the paper (added after Molin's remark of 10:30)

The comparison figures of Section 2.3 (`claude_small_patch/results/`) now show only the new runs - star re-run
(solid, filled circles) and single-column patch (dashed, open squares) - with no figure titles (the panel labels
'star patch' / 'single-column patch' and '(a)'-'(d)' remain); the paper's own data enter the tables only. All
figures of this study and the paper's figures below share one style (`claude_small_patch/paper_style.py`: viridis
ramp cut at 0.82 so that the lightest member prints, dark = first member of a sweep; 10 pt fonts, 7 pt legends
where crowded; light grid; sizes chosen so that the fonts print close to nominal in the SIAM one-column layout:
`SINGLE` = 4.6 x 3.1 in for a figure at 0.8-0.85 linewidth, `HALF` = 2.8 x 2.4 in for the 0.48-linewidth nonlinear
panels, `WIDE` = 5.8 x 2.9 in for the two-panel history figures, `QUAD` = 6.4 x 5.2 in for the four-panel figures).

The paper's ten iteration figures were regenerated in that style with the same file names, so they are drop-in
replacements (`restyle_paper_figures.py`):

| folder under `../Patch-Preconditioner-for-Boussinesq/Image_restyled/` | data | use |
|---|---|---|
| `.` (root) | the paper's own folders (`data_linear_*_deg2`, `data_nonlinear_*_deg2_ew`): the content of figs. 4-12 unchanged, style only; plus `column_patch_sweeps.png` (fig. 14, same data as `Image_patch/`) | replace `Image_new/`, `Image_nonlinear/`, `Image_patch/` in the tex by `Image_restyled/` (one `sed`), or copy the files over |
| `star_rerun/` | the 2026-10-16 star re-run of the same sweeps (identical curves except near the iteration cap) | not needed for the paper; shows that the code at HEAD reproduces the figures |
| `single_column/` | the single-column patch, same sweeps, same figure names | the column-patch versions of figs. 4-12, for the single-column section |
| `corrected/` under each of the above | the same figures with the iteration counts as counts: the linear shift sweep plots n (the paper's plots n + 1), the nonlinear shift sweep plots the total KSP count over the Newton steps (the paper's plots total + number of Newton steps), and the linear histories keep their last point (the converged state, which the paper's figures drop) | if the text is changed to the corrected conventions |

Conventions of the drop-in set are those of the paper's figures, figure by figure (linear histories drop the final
point; the nonlinear histories keep every Newton iterate including the converged one, except
`snes_error_AR_dt5000.0`, whose original drops it; the two shift-sweep figures use the old counting of
`plot_shift.py`), so that nothing but the style changes. Note that `plot_dx.py`, `plot_dz.py`, `plot_ar.py` and
`plot_shift.py` at HEAD all drop the converged SNES point (commit 82cd238 of 2026-07-08, after the paper's figures
were made), so Molin's own re-plot of the paper's data would not reproduce four of the five nonlinear figures either;
`corrected/` uses that one convention throughout. Legends read $N_x$, $N_z$, $L_z$ (with $\varepsilon$ in the linear
aspect-ratio figure) and $C_1$ instead of `nx=`, `AR=` and the mislabelled `K=`. Each restyled figure was compared
with its original by an independent reader (curve by curve, end points, legend entries); Section 7.1.

To switch the paper to the restyled set:

```bash
cd ../Patch-Preconditioner-for-Boussinesq
sed -i 's#Image_new/#Image_restyled/#g; s#Image_nonlinear/#Image_restyled/#g; s#Image_patch/#Image_restyled/#g' ex_article.tex
```

(the `\includegraphics` widths can stay; the half-width panels were drawn at 2.8 in for the 0.48-linewidth
minipages, the single figures at 4.6 in for 0.8-0.85 linewidth, fig. 14 at 6.4 in for the full linewidth.)

### 7.1 Independent check of the restyled figures

Eleven agents compared each restyled figure with its original pixel by pixel (curve counts, points per curve,
end points read off the tick calibration, axis ranges, legend entries) and judged legibility at the printed width;
a twelfth checked the comparison figures of `results/`. Findings and what was done about them:

* Ten of the eleven figures showed the same data within reading accuracy (1-3 %, one pixel on the log axes); the
  one difference was `snes_error_AR_dt5000.0`: the paper's figure has three points per curve (it drops the
  converged iterate), the first restyled version four. Fixed (drop_last for that figure).
* Legends on data: `snes_error_shift_dt500.0` (seven entries at half width, first placed lower left over the fast
  curves) and panel (a) of fig. 14 (the two-column legend hid the maximum at dt = 5e4-1e5 s) and its panel (c)
  (the K = 100 peak). Fixed (upper right in the empty corner; y-ranges extended to make room; short labels).
* Print size: the first sizes (6.0 / 3.4 / 7.6 in) would have printed the 8 pt legends at 4.5-6 pt. Fixed by the
  print-true sizes above.
* The lightest viridis colour had contrast 1.6:1 on white; the ramp now stops at 0.82.
* In the two-panel history figures the single-column panel first used the star style; it now uses dashed lines and
  open squares.
* Inherent in the data, same in the originals, left as is: curves that coincide (N_z = 60 under N_z = 80 in
  `snes_error_dz_dt50.0`, N_x = 40/60 under N_x = 80 in `error_dx_Robust_dt50.0`), and the last labelled x tick at
  10^4 s while the nonlinear shift sweep extends to 5e4 s.

Not re-checked by an agent after these changes; I looked at every changed figure myself.

---

The follow-up of 2026-10-19 (the paper edits and the Table-5 study) is recorded in `claude_paper_update_and_table5.md`.
