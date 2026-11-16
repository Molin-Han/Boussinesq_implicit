# Single-column (Vanka) patch versus star patch

Study notes for proofreading. Written 2026-09-29.
Everything referred to below lives in `claude_column_patch/` next to this file, and in
`../data_claude_paper_nonlinear_{dx,dz,ar,dt}_patch{star,vanka}/`.

Section 8 lists what is solid and what is not. Read it before relying on any statement here.

---

## 0. Answers in short

**1. Is the idea correct mathematically?** Partly.

* The *premise* is not correct. The discrete divergence does not degrade when the domain gets thinner:
  its matrix, and therefore its kernel, is the same for every aspect ratio (Section 2.2).
  The horizontal part of the grad-div term is negligible only when
  $\gamma = c/\Delta x^2 = \Delta t/(2\delta\Delta x^2) \lesssim 1$, which for $\delta = K\delta_{thr}$ means
  $K \gtrsim N_x^2\max(1, 4\varepsilon^2 m)/(4\pi^2)$, about 160 on the nominal mesh. This is a condition on the shift
  and the horizontal resolution, not on thinness.
* The *conclusion* (a single-column patch can replace the star patch) is nevertheless valid for degree 2, for a
  different reason: the quadratic bubbles in $x$ span exactly the part of the divergence-free space that the coarse
  level misses, and they live in single columns (Section 2.4). For degree 1 it fails.
* The price is a larger constant: between 2 and 12.3 in the two-level analysis, against about 1.3 for the star patch.
  It is *largest* in the thin, weakly stratified limit, and there it also grows with the number of levels (Section 2.5).
  So thin domains are the least favourable case for the column patch, not the most favourable.

**2. Code.** `--patch vanka` puts `firedrake.ASMVankaPC` with `pc_vanka_construct_dim 1` inside `mg_levels`.
The patch is what was intended (one cell column plus the vertical facets on both sides, 364 unknowns against 608),
serial and parallel. `--patch star` (default) reproduces the committed solver exactly (Section 3).

**3. Experiments** (440 sweep runs, 60 timing runs, Sections 4 and 5).

| | star | single column |
|---|---|---|
| Converges wherever the other does | yes | yes (identical set of failures, all at $K \ge 100$) |
| Total FGMRES iterations, $K = 1$ | 3 to 40 | 10 to 50 |
| Ratio to star, $K = 1$: median / max | | dx 1.39 / 5.0, dz 1.63 / 6.5, ar 1.67 / 10.7 |
| Ratio to star, $K \ge 100$ | | 0.94 to 1.06 (identical counts in 6 of 10 runs) |
| Multigrid alone, random right-hand side | 2 | 6 to 20 (bounded in every parameter tested) |
| Time per FGMRES iteration (nominal mesh) | 0.28 to 0.31 s | 0.20 to 0.22 s (about 30 % less) |
| Time of one time step | | 0.73 to 0.93 of star where counts are close; 1.2 to 4.1 of star at small $\Delta t$, thin domain, fine $\Delta z$ |
| Iteration counts reproducible run to run | yes | no, at $\Delta t \le 100$ s (e.g. 13 to 20 for the same command) |

So the parameter robustness survives in the sense that the counts stay bounded and nothing new fails, but the solve is
faster only where the iteration counts of the two patches are within about 1.4 of each other, that is for
$\rho_1 = mL_z^2/(\pi^2\Delta x^2) \gtrsim 10$ (large $\Delta t$) or for $K \gtrsim 100$.

---

## 1. What I changed and what I created

### 1.1 Edited files (two)

**`Nonlinear_Boussinesq_Irk_SC.py`** (tracked; `git diff` shows the change, nothing is committed):

```diff
+parser.add_argument('--patch', type=str, default='star', choices=['star', 'vanka'], help='...')
-solver_name = 'MG ASMStar'
+solver_name = 'MG ASMStar' if args.patch == 'star' else 'MG ASMVanka'
+# ! Patch for the ASM smoother on the MG levels.
+# ! The topology_dm of the extruded mesh is the base interval mesh, so every patch is made of whole columns.
+if args.patch == 'star':
+    patch_params = { "pc_type": "python", "pc_python_type": "firedrake.ASMStarPC", "pc_star_construct_dim": 0,
+                     "pc_star_sub_sub_pc_type": "lu", 'pc_star_sub_sub_pc_factor_mat_ordering_type': 'rcm',
+                     'pc_star_sub_sub_pc_factor_reuse_ordering': None }
+else:
+    patch_params = { "pc_type": "python", "pc_python_type": "firedrake.ASMVankaPC", "pc_vanka_construct_dim": 1,
+                     "pc_vanka_sub_sub_pc_type": "lu", 'pc_vanka_sub_sub_pc_factor_mat_ordering_type': 'rcm',
+                     'pc_vanka_sub_sub_pc_factor_reuse_ordering': None }
 'mg_levels': { 'ksp_type': 'gmres', 'ksp_max_it': 6,
-    "pc_type": "python", "pc_python_type": "firedrake.ASMStarPC", ... (the six star lines)
+    **patch_params,
 },
```

**`Bash/bash_claude_paper.py`** (untracked, from the 2026-09-20 session): new option `--patch {star,vanka}`.
It appends `--patch vanka` to the nonlinear commands, appends nothing for `star` (so the earlier commands are
reproduced character for character), and refuses `--mode linear --patch vanka` because `LB_time_slice.py` has no such flag.

### 1.2 New files

| Path | What |
|---|---|
| `claude_column_patch_study.md` | this file |
| `claude_column_patch/sweeps/` | `compare_patches.py`, `quick_table.py`, `compare_iterations.md`, `compare_{dx,dz,ar,dt}.png`, `collapse_rho1.png` |
| `claude_column_patch/timing/` | `run_timing.sh`, `parse_timing.py`, `runs/` (60 runs, `timing_summary.md`, `timing_all_runs.csv`) |
| `claude_column_patch/my_checks/` | `check_two_level_constant.py` (my independent check of Section 2.5), `degree1_runs/` |
| `claude_column_patch/theory_kernel/` | model of the discretisation through vertical modes, kernel counts, constants (by an agent) |
| `claude_column_patch/toy_model/` | lowest-order model problem, 5250-configuration spectra, multigrid sweeps (by an agent) |
| `claude_column_patch/firedrake_algebraic/` | operators and kernels of the real discretisation (by an agent) |
| `claude_column_patch/firedrake_inner_solve/` | `inner_solve.py`: the field-1 multigrid alone; `results_sweep_*.csv` |
| `claude_column_patch/code_verification/` | option comparison, patch enumeration, serial/parallel runs (by an agent) |
| `claude_column_patch/theory_asymptotic/` | scripts only, no results written (agent cut off) |
| `test_claude_paper_nonlinear_{dx,dz,ar,dt}_patch{star,vanka}.sh` | written by `run_claude_paper.sh`, serial form of the sweeps |
| `../data_claude_paper_nonlinear_{dx,dz,ar,dt}_patch{star,vanka}/` | the sweep data, 2.2 GB in total (mostly VTK output inside `jobs/`) |

Memory (in `~/.claude/projects/-home-mh2520-molin-phd-yr2-Boussinesq-implicit/memory/`): your name and the
short-answer preference, plus the four project notes copied from the parent folder's memory.

### 1.3 Not touched

`utils.py`, `LB_time_slice.py`, `run_claude_paper.sh`, `plot_claude_paper.py`, the earlier data folders. No commit, no push.

---

## 2. Mathematics

### 2.1 Setting

Field 1 of the auxiliary operator, linearised about rest with $f = 0$: $u_y$ decouples (mass matrix), and eliminating
$b$ leaves

$$a(u, w) = (u_x, w_x) + m\,(u_z, w_z) + c\,(\nabla\!\cdot u, \nabla\!\cdot w),\qquad
m = 1 + \frac{\Delta t^2N^2}{4},\qquad c = \frac{\Delta t}{2\delta} = \frac{C_M^2}{K}.$$

Checked against the running solver on the real discretisation
(`firedrake_algebraic/results_01_verify_*.txt`): the matrix extracted from `fieldsplit_1` equals the hand-written
form to $4\times10^{-15}$, the Schur complement after eliminating $b$ equals the form above to $3\times10^{-16}$,
the $u$–$u_y$ and $u_y$–$b$ blocks are exactly zero.

$V_h$ is the Raviart–Thomas space of degree $k$ on rectangles, periodic in $x$, $u\cdot n = 0$ on top and bottom.
The multigrid coarsens in $x$ only by factors of 2; every level has the same $N_z$. Both patches consist of whole
columns, because the DMPlex of an extruded mesh is the base mesh.

### 2.2 The premise: does the divergence degrade in a thin domain?

No, in two senses.

**Algebraically.** In the basis Firedrake uses (Piola-mapped), the divergence matrix does not depend on the geometry:

| compared with $L_x = 300$ km, $L_z = 4$ km | $\max\lvert B - B_{ref}\rvert$ | $\max\lvert B\rvert$ |
|---|---|---|
| $L_z = 250$ m | $4\times10^{-26}$ | 0.577 |
| $L_z = 12$ km | $3\times10^{-16}$ | 0.577 |
| $L_x = L_z = 1$ | $2\times10^{-16}$ | 0.577 |

(`firedrake_algebraic/results_02_kernel_capture.txt`, degree 2, 16 columns $\times$ 6 layers.) So
$N_h = \ker(\nabla\!\cdot)$ has the same dimension, $k N_x(kN_z - 1) + 1$, and the same elements for every aspect
ratio. Thinness changes only the mass matrix, i.e. the inner product in which the kernel is measured.

**In magnitude.** For mesh-scale fields the horizontal part of $c\,(\nabla\!\cdot u,\nabla\!\cdot w)$ compares with the
mass term of $u_x$ as $\gamma = c/\Delta x^2$. With $\delta = K\delta_{thr}$,

$$\gamma = \frac{c}{\Delta x^2} = \frac{N_x^2\,\max(1,\,4\varepsilon^2 m)}{4\pi^2 K},$$

which does not contain $\varepsilon$ at all in the weakly stratified regime. On the nominal mesh ($N_x = 80$, $K = 1$)
$\gamma = 162$: the horizontal part is far from negligible. It becomes negligible for $K \gtrsim 160$.

The same file records what would happen if the divergence really were $\partial_z u_z$ alone: then both patches
capture the whole kernel (span of patch kernels = 128 = dimension of the kernel). So the reasoning from the premise
to the conclusion is sound; it is the premise that does not hold.

### 2.3 What each patch contains of the kernel

$N_h = \mathrm{curl}\,\Psi_h^0 \oplus \mathrm{span}\{(1,0)\}$ with $\Psi_h = CG_k(x)\otimes CG_k(z)$ and $\Psi_h^0$ its
members vanishing on top and bottom. On $N_h$ the form is an anisotropic Laplacian of the stream function,
$a(\mathrm{curl}\,\psi, \mathrm{curl}\,\chi) = (\partial_z\psi, \partial_z\chi) + m(\partial_x\psi, \partial_x\chi)$.

| | star of base vertex $i$ | closure of base cell $j$ (single column) |
|---|---|---|
| stream functions with curl inside the patch space | $\varphi(x)g(z)$, $\varphi\in CG_k(x)$ supported in the two cells at $x_i$ | $\beta(x)g(z)$, $\beta$ a bubble of cell $j$ |
| dimension, general | $(2k-1)(kN_z - 1)$ | $(k-1)(kN_z - 1)$ |
| degree 2, $N_z = 20$ | 117 | 39 |
| degree 1 | $N_z - 1$ | **0** |
| sum over all patches | all of $N_h$ except the uniform flow | **half** of it (degree 2), **nothing** (degree 1) |

Why the column patch has so little: a facet basis function lives in both neighbouring cells, so a field of the patch
space with non-zero flux through a lateral facet has non-zero divergence in the neighbouring column, where the patch
has nothing to compensate it. Only fields with zero flux through both lateral facets can be divergence free.

Measured on the real discretisation (`firedrake_algebraic/results_02_kernel_capture.txt`), degree 2, $16\times6$:
kernel 353; star patches span 352; column patches span 176. Degree 1, $16\times6$: kernel 81; star 80; column 0.
All counts match the formulas.

### 2.4 Why the column patch nevertheless works for degree 2

The coarse level must be counted. For quadratics in $x$ and coarsening by 2,

$$CG_2(x; h) = CG_2(x; 2h) \oplus \mathrm{span}\{\text{cell bubbles of the fine mesh}\},$$

because the fine vertices are exactly the nodes of $CG_2(2h)$ (coarse vertices and coarse midpoints): subtract from
$U$ the coarse quadratic $I_HU$ that interpolates it at the fine vertices, and what remains vanishes at every fine
vertex, i.e. is a sum of fine bubbles. Tensorised with $CG_2(z)$:

$$N_h = N_{2h} \oplus \sum_j \big(N_h\cap V_j^{col}\big)\qquad (k = 2).$$

This is the kernel condition of the subspace-correction theory (Schöberl; Lee–Wu–Xu–Zikatanov), satisfied by
*coarse level + column patches*. The same argument works for every $k\ge2$. For $k = 1$ the complement of
$CG_1(2h)$ in $CG_1(h)$ is spanned by the hats of the odd vertices, which live on two cells: no single-column
function exists, and the condition fails.

Measured with Firedrake's own prolongation (`theory_kernel/02_firedrake_crosscheck.txt`, $8\times3$):

| degree | $\dim N_h$ | $\dim(N_{2h} + \text{column patch kernels})$ | $\dim(N_{2h} + \text{star patch kernels})$ |
|---|---|---|---|
| 1 | 17 | **9** | 17 |
| 2 | 81 | **81** | 81 |

So your instinct that a smaller patch can carry a local basis is right for the element you use, but the basis is the
set of $x$-bubbles, and it has nothing to do with the domain being thin.

### 2.5 The price: the constant of the splitting

Because patches are whole columns, the 2D operator separates exactly into vertical modes $\sigma_n$
(checked: the union of the mode spectra equals the 2D spectrum to $10^{-14}$,
`theory_kernel/03_mode_decoupling_check.txt`). In mode $\sigma$ the divergence-free fields are parametrised by
$U\in CG_k(x)$ with energy

$$S_\rho(U) = \lVert U\rVert^2 + \rho\,h^2\lVert U'\rVert^2,\qquad
\rho = \frac{m}{(\sigma\Delta x)^2},\qquad \rho_1 = \frac{mL_z^2}{\pi^2\Delta x^2}\ \text{ for the gravest mode.}$$

The two-level constant of the column patch for $c\to\infty$ is
$C(\rho) = \sup_U\,[S_\rho(I_HU) + S_\rho(U - I_HU)]/S_\rho(U)$, and $\lambda_{min}\to 1/C(\rho)$.

| $\rho$ | 0 | 0.001 | 0.01 | 0.12 | 1 | 10 | $\infty$ |
|---|---|---|---|---|---|---|---|
| $C(\rho)$ | 12.28 | 11.71 | 8.50 | 3.35 | 2.19 | 2.02 | 2.00 |

I recomputed this table with a script written independently of the agent's code
(`my_checks/check_two_level_constant.py`); the two agree to all printed digits.

Reading: $\rho\to0$ means the energy is $L^2$-like in $x$, which is the thin domain with weak stratification
($L_z\sqrt m \ll \Delta x$); $\rho\to\infty$ means $H^1$-like in $x$ (thick domain or large $\Delta t$). Interpolation
at vertices is a poor operator in $L^2$ and a good one in $H^1$, hence 12.3 against 2.

Several levels turn this into the hierarchical-basis splitting in $x$, which is stable in $H^1$ and not in $L^2$:

| levels | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| $\rho = 0$ | 12.3 | 15.7 | 36.0 | 73.6 |
| $\rho = 1$ | 2.19 | 3.16 | 4.91 | 9.07 |
| $\rho = \infty$ | 2.00 | 2.46 | 2.67 | 2.78 |

(my script; the agent's multilevel $\lambda_{min}$ at $\rho\approx0$, 0.083, 0.064, 0.028, 0.014, are the reciprocals
to within 2 %.) Your code has 3 levels. With more levels the column patch would lose robustness in thin domains.

For comparison, the star patch: two-level $\lambda_{min}$ between 0.79 and 0.83 for degree 2 in every sweep of
`theory_kernel/04_spectra_parameter_study.txt`, independent of $\gamma$, $L_z$, $N_x$, $N_z$, and of the number of
levels.

Additive two- and three-level $\lambda_{min}$ at the nominal configuration, degree 2 (same file):

| $K$ | $\gamma$ | star 3 levels | column 1 level | column 2 levels | column 3 levels |
|---|---|---|---|---|---|
| 0.01 | 16211 | 0.829 | $4.1\times10^{-6}$ | 0.0816 | 0.0639 |
| 1 | 162 | 0.829 | $4.1\times10^{-4}$ | 0.0820 | 0.0657 |
| 100 | 1.6 | 0.829 | 0.039 | 0.118 | 0.137 |
| 1000 | 0.16 | 0.829 | 0.275 | 0.331 | 0.341 |

The column patch alone is not robust ($\lambda_{min}\propto1/\gamma$); with the coarse level it is, at a constant
about 12 times worse than the star.

### 2.6 Degree 1

Lowest-order model, 5250 configurations (`toy_model/05_analysis_summary.txt`): the condition number of the
column-patch preconditioner follows $\mathrm{cond}\approx\gamma\,\max(1.5,\ 13.5/(1+\eta))$ with
$\eta = m\Delta z^2/\Delta x^2$, to within 13 % for $\gamma>10$ (fitted slope in $\gamma$: 0.998). About 45 % of all
eigenvalues are below a tenth of the largest, against 0.5 % for the star. Making the domain thinner never improved
it in the weak regime (0 of 3300 halving steps), and did in the strong regime (99 % of 900 steps).

### 2.7 Predictions, and what the sweeps show

| Prediction | Observed |
|---|---|
| Counts equal to star when $\gamma\lesssim1$, i.e. $K\gtrsim160$ | $K = 100$ and $1000$: ratio between 0.94 and 1.06 in all 10 runs where both converge, identical in 6; $K = 10$: 1.0 to 1.67 (Section 4.2) |
| Bounded as $K\to0$ for degree 2 | multigrid alone: yes, 17 and 20 at $K = 1, 0.1$ (Section 4.4). **Full solver: not clear.** At $K = 0.1$ the column counts are 12 to 55 (star 4 to 16), more than at $K = 1$ (13 to 31); see the note after the dt table |
| Not bounded for degree 1 | multigrid alone: 19 at $K = 1$, no convergence in 200 at $K = 0.1$; full solver at $K = 1$: 56, 78, 99 for $N_x = 10, 20, 40$ (star 11, 13, 13) |
| Penalty governed by $\rho_1$, largest for $\rho_1\to0$ | median ratio 5.2, 3.0, 2.4, 1.5, 1.16, 1.11 in the six decades of $\rho_1$ (Section 4.3) |
| Thinner domain: worse | $\Delta t = 50$: 17, 18, 24, 25, 29, 29, 32 for $L_z = 12 \dots 0.25$ km, star 13 down to 3 |
| Finer $\Delta x$: better ($\rho_1$ grows) | $\Delta t = 50$: 23, 24, 19, 16, 15 for $N_x = 10\dots80$ |
| Large $\Delta t$: close to star ($m$, hence $\rho_1$, large) | ratio $\le1.2$ for $\Delta t\ge5000$ on the nominal geometry |
| Independent of $N_z$ | **not confirmed in the full solver**: 20, 26, 27, 32, 39 at $\Delta t = 50$; confirmed for the multigrid alone (16, 17, 17, 17) |

---

## 3. Code-level verification

**The flag does what it says.** `code_verification/01_static_options_check.txt`: the option set of the committed
script and of the working tree with `--patch star` are identical (41 options, none differing), also with `--timing`
(39). `--patch vanka` differs in exactly five entries: `pc_python_type` and the four `pc_star_*` keys replaced by
their `pc_vanka_*` counterparts. The solver views (`ksp_view_{star,vanka}_nominal_2ranks.txt`) show
`Python: firedrake.ASMVankaPC` under `mg_levels_1_` and `mg_levels_2_`, patch matrices of 364 rows.

**Star path unchanged.** Residual histories of `--patch star` and of the committed script, serial: identical in every
printed digit. The four new star sweeps reproduce the sweeps of 2026-09-21 in all 220 runs they have in common
(total iterations and status).

**Patch composition** (degree 2; enumerated with Firedrake's own `build_*_indices` and compared with the patches
built from the definition; `code_verification/runs/patch_dofs_*/report.txt`):

| field, entity | single column | star |
|---|---|---|
| $u$, vertical facets | $4N_z$ (both sides) | $2N_z$ (the middle ones) |
| $u$, horizontal facets | $2(N_z+1)$ | $4(N_z+1)$ |
| $u$, cell interiors | $4N_z$ | $8N_z$ |
| $u_y$ | $4N_z$ | $8N_z$ |
| $b$ | $4N_z+2$ | $8N_z+4$ |
| total | $18N_z + 4$ = 364 | $30N_z + 8$ = 608 |
| LU factor non-zeros ($N_z = 20$) | 21 060 | 53 468 |

Every column patch runs from bottom to top, contains one cell column, and the vertical facets at exactly its two
base vertices. Each vertical-facet unknown is in 2 column patches, everything else in 1 (star: the other way round).

**Parallel.** The views report 44 / 84 blocks for the column patch on levels with 40 / 80 cells (star 46 / 86). The
extra ones are *empty* index sets, one per ghost entity: Firedrake loops over all local entities and skips the ghost
seeds. With 2 and 4 ranks: 40 / 80 non-empty patches, all identical to the definition, none missing, none
duplicated, none incomplete (`runs/patch_dofs_np{2,4}_*`). The warning
`Mesh overlap depth of 1 does not support vanka-patches` is generic (it assumes vertex-centred Vanka patches) and
does not apply to `construct_dim 1` here: the closure of an owned cell only needs its two vertices, which the
overlap of depth 1 provides.

**A trap worth knowing.** PETSc options given on the command line do *not* override entries of the
`solver_parameters` dictionary here, because the solver has an automatically generated prefix (`firedrake_1_`)
and petsctools then ignores the command line for keys that the dictionary sets. My first attempt,
`-firedrake_1_shiftedschurpc_fieldsplit_1_mg_levels_pc_python_type firedrake.ASMVankaPC`, silently ran the star patch.
Options that are *not* in the dictionary (`-firedrake_1_snes_atol`, `-firedrake_1_snes_ksp_ew`) do work. Hence the flag
in the script.

---

## 4. Numerical experiments: iteration counts

### 4.1 Set-up

The four sweeps of `claude_paper_tests_README.md`, run for both patches with the same commands except `--patch vanka`:
nonlinear code, $f = 0$, no Eisenstat–Walker, `rtol` $10^{-6}$, `ksp_atol` $10^{-9}$, `snes_atol` $10^{-9}$,
`maxit` 100, 2 MPI ranks, $\delta = K\delta_{thr}$, second time step monitored. $K\in\{0.1, 1, 10, 100, 1000\}$ in the
dt sweep, $K = 1$ elsewhere. 220 runs per patch. Numbers are total FGMRES iterations over all Newton steps,
written `star / column`.

### 4.2 Tables

| sweep | runs compared | median ratio | mean | max | not converged, star | not converged, column |
|---|---|---|---|---|---|---|
| dx | 50 | 1.39 | 1.72 | 5.00 | 0 | 0 |
| dz | 50 | 1.63 | 2.18 | 6.50 | 0 | 0 |
| ar | 70 | 1.67 | 2.55 | 10.67 | 0 | 0 |
| dt | 40 | 1.11 | 1.76 | 5.50 | 10 of 50 | 10 of 50 (the same ten) |

**dx** ($N_x$ of the coarse mesh; 4 times as many columns on the finest level)

| $\Delta t$ | 10 | 20 | 40 | 60 | 80 |
|---|---|---|---|---|---|
| 1 | 4 / 19 | 4 / 20 | 4 / 13 | 4 / 13 | 4 / 10 |
| 10 | 6 / 18 | 6 / 18 | 6 / 14 | 6 / 16 | 6 / 14 |
| 50 | 13 / 23 | 13 / 24 | 11 / 19 | 11 / 16 | 11 / 15 |
| 100 | 12 / 25 | 12 / 24 | 10 / 20 | 10 / 18 | 13 / 14 |
| 500 | 16 / 25 | 17 / 29 | 16 / 27 | 16 / 25 | 16 / 22 |
| 1000 | 20 / 29 | 20 / 29 | 20 / 28 | 21 / 27 | 20 / 25 |
| 5000 | 28 / 31 | 29 / 31 | 29 / 30 | 27 / 29 | 27 / 29 |
| 10000 | 28 / 30 | 28 / 29 | 27 / 28 | 25 / 27 | 25 / 27 |
| 50000 | 19 / 21 | 19 / 22 | 19 / 23 | 19 / 23 | 19 / 23 |
| 100000 | 19 / 21 | 18 / 20 | 18 / 21 | 18 / 19 | 18 / 19 |

**dz** ($N_z$)

| $\Delta t$ | 10 | 20 | 40 | 60 | 80 |
|---|---|---|---|---|---|
| 1 | 6 / 17 | 4 / 18 | 4 / 21 | 4 / 21 | 4 / 26 |
| 10 | 6 / 18 | 6 / 20 | 6 / 18 | 6 / 22 | 6 / 31 |
| 50 | 11 / 20 | 13 / 26 | 11 / 27 | 11 / 32 | 11 / 39 |
| 100 | 12 / 20 | 12 / 24 | 13 / 29 | 13 / 50 | 13 / 40 |
| 500 | 16 / 25 | 17 / 29 | 16 / 38 | 16 / 42 | 16 / 41 |
| 1000 | 20 / 26 | 20 / 29 | 20 / 32 | 21 / 32 | 20 / 36 |
| 5000 | 28 / 30 | 29 / 31 | 29 / 32 | 29 / 33 | 31 / 36 |
| 10000 | 28 / 28 | 28 / 29 | 27 / 30 | 29 / 33 | 29 / 33 |
| 50000 | 19 / 23 | 19 / 22 | 19 / 22 | 19 / 23 | 19 / 23 |
| 100000 | 18 / 20 | 18 / 20 | 18 / 21 | 18 / 21 | 18 / 22 |

**ar** ($L_z$ in m)

| $\Delta t$ | 12000 | 8000 | 4000 | 2000 | 1000 | 500 | 250 |
|---|---|---|---|---|---|---|---|
| 1 | 8 / 11 | 6 / 13 | 4 / 18 | 4 / 18 | 5 / 28 | 5 / 27 | 3 / 27 |
| 10 | 8 / 13 | 9 / 14 | 6 / 18 | 4 / 29 | 6 / 27 | 5 / 26 | 3 / 20 |
| 50 | 13 / 17 | 10 / 18 | 13 / 24 | 10 / 25 | 6 / 29 | 6 / 29 | 3 / 32 |
| 100 | 13 / 20 | 12 / 21 | 12 / 24 | 11 / 27 | 7 / 33 | 6 / 31 | 5 / 27 |
| 500 | 18 / 25 | 18 / 25 | 17 / 29 | 14 / 33 | 13 / 40 | 11 / 35 | 11 / 39 |
| 1000 | 23 / 27 | 21 / 28 | 20 / 29 | 17 / 33 | 15 / 38 | 13 / 36 | 11 / 35 |
| 5000 | 23 / 24 | 26 / 28 | 29 / 31 | 27 / 35 | 27 / 35 | 20 / 35 | 18 / 37 |
| 10000 | 22 / 23 | 23 / 24 | 28 / 29 | 32 / 34 | 27 / 38 | 27 / 39 | 22 / 38 |
| 50000 | 15 / 19 | 18 / 19 | 19 / 22 | 26 / 28 | 32 / 33 | 37 / 39 | 34 / 44 |
| 100000 | 15 / 19 | 15 / 19 | 18 / 20 | 19 / 23 | 26 / 28 | 32 / 33 | 40 / 44 |

In thin domains at small $\Delta t$ the star patch needs 3 to 6 iterations and the column patch 20 to 32.

**dt** ($K$; `maxit` = a linear solve reached 100 iterations, `div` = Newton diverged)

| $\Delta t$ | 0.1 | 1 | 10 | 100 | 1000 |
|---|---|---|---|---|---|
| 1 | 4 / 12 | 4 / 13 | 6 / 10 | 10 / 10 | 18 / 18 |
| 10 | 4 / 17 | 6 / 18 | 9 / 10 | 16 / 16 | 51 / 51 |
| 50 | 8 / 44 | 13 / 24 | 18 / 21 | 62 / 62 | maxit / maxit |
| 100 | 9 / 30 | 12 / 23 | 21 / 23 | 65 / 61 | maxit / maxit |
| 500 | 11 / 55 | 17 / 29 | 32 / 34 | 106 / 101 | maxit / maxit |
| 1000 | 12 / 53 | 20 / 29 | 40 / 43 | 125 / 132 | maxit / maxit |
| 5000 | 16 / 26 | 29 / 31 | 58 / 58 | maxit / maxit | maxit / maxit |
| 10000 | 16 / 23 | 28 / 29 | 60 / 61 | maxit / maxit | div / div |
| 50000 | 12 / 20 | 19 / 22 | 32 / 32 | 82 / 84 | maxit / maxit |
| 100000 | 12 / 18 | 18 / 20 | 27 / 27 | 55 / 55 | maxit / maxit |

Note on small shifts. Going from $K = 1$ to $K = 0.1$ helps the star patch (13 to 8 at $\Delta t = 50$) and hurts the
column patch (24 to 44), almost entirely in the first Newton step (35, 36, 34 iterations at $\Delta t = 50, 500, 1000$).
A plausible reason, which I did not verify: the auxiliary operator recovers the pressure as $p - \delta^{-1}\nabla\!\cdot u$,
so an error of the field-1 solve is amplified by $1/\delta$; the star multigrid is accurate enough for this not to
matter, the column multigrid is not. Two single runs made by an agent at $K = 0.1$ and $0.01$ gave 23 and 25, which
I could not reconcile with the 44 of the sweep; see Section 4.5 on reproducibility.

Figures: `claude_column_patch/sweeps/compare_{dx,dz,ar,dt}.png`.

### 4.3 Collapse against $\rho_1$

All 170 runs at $K = 1$ (dx, dz, ar), binned by $\rho_1 = mL_z^2/(\pi^2\Delta x^2)$ with $\Delta x$ of the finest level:

| $\rho_1$ | runs | median column/star | max | median extra iterations | max extra |
|---|---|---|---|---|---|
| below 0.01 | 13 | 5.20 | 10.67 | 23 | 29 |
| 0.01 to 0.1 | 12 | 3.04 | 7.25 | 16 | 27 |
| 0.1 to 1 | 49 | 2.36 | 6.50 | 12 | 37 |
| 1 to 10 | 27 | 1.52 | 3.25 | 8 | 16 |
| 10 to 100 | 19 | 1.16 | 1.41 | 4 | 11 |
| above 100 | 50 | 1.11 | 1.27 | 2 | 4 |

Figure: `claude_column_patch/sweeps/collapse_rho1.png`. The scatter inside a bin is mostly the dz sweep (next item).

### 4.4 The multigrid of field 1 alone

`firedrake_inner_solve/inner_solve.py` solves only the field-1 problem with FGMRES to $10^{-8}$, preconditioned by
exactly the multigrid of the solver script. It was checked against the multigrid object of the running solver: same
iteration count (17) and residual history to $10^{-4}$ (`verification_vanka.json`). Serial, degree 2, nominal
configuration unless varied. `random` = random right-hand side, `bump` = localised smooth one.

| varied | values | star, random | column, random | column, bump |
|---|---|---|---|---|
| $K$ | 0.1, 1, 10, 100, 1000 | 2, 2, 2, 2, 1 | 20, 17, 8, 5, 3 | 8, 7, 4, 2, 2 |
| $\Delta t$ | 1, 50, 1000, $10^4$, $10^5$ | 2 each | 17, 17, 16, 9, 6 | 7, 7, 4, 2, 3 |
| $N_x$ (coarse) | 10, 20, 40, 80 | 2 each | 11, 17, 19, 19 | 7, 7, 5, 4 |
| $N_z$ | 10, 20, 40, 80 | 2 each | 16, 17, 17, 17 | 6, 7, 7, 8 |
| $L_z$ (m) | 12000, 4000, 1000, 250 | 2 each | 16, 17, 17, 18 | 5, 7, 12, 15 |

The star multigrid is close to an exact solver (one application reduces a smooth residual by about $10^{-11}$).
The column multigrid is bounded in every parameter but needs about eight times as many iterations, and for the
localised right-hand side it degrades as the domain thins (5 to 15).

**More smoothing does not rescue it** (column patch, nominal, random):

| GMRES steps in the smoother | 6 | 12 | 20 |
|---|---|---|---|
| iterations | 17 | 9 | 8 |
| patch solves in total | 38 080 | 37 440 | 53 760 |

**Degree 1** (real discretisation, nominal, random right-hand side): star 3, 3, 2, 2 for $K = 0.1, 1, 10, 100$;
column: no convergence in 200, 19, 13, 6.

Degree 1 in the full solver ($K = 1$, $\Delta t = 50$, 2 ranks, `my_checks/degree1_runs/`):

| $N_x$ (coarse) | 10 | 20 | 40 |
|---|---|---|---|
| star | 11 | 13 | 13 |
| column | 56 | 78 | 99 |

### 4.5 Reproducibility

The nominal configuration occurs in all four sweeps with the same command. Star: identical counts in all four, at
every $\Delta t$. Column: identical for $\Delta t\ge500$, but

| $\Delta t$ | 1 | 10 | 50 | 100 |
|---|---|---|---|---|
| column, the four sweeps | 20, 18, 18, 13 | 18, 20, 18, 18 | 24, 26, 24, 24 | 24, 24, 24, 23 |

and three identical timing runs at $N_z = 80$ gave 46, 35, 38; at $L_z = 250$ m 27, 30, 27. Serial and 2-rank runs
differ in the same way (25 against 24 to 26). The patches are the same (Section 3), so this is sensitivity of a slowly
converging iteration with a Krylov method inside the preconditioner to rounding-level differences. Treat individual
column counts at small $\Delta t$ as $\pm20\,\%$.

---

## 5. Timing

`claude_column_patch/timing/run_timing.sh`: `--timing`, 2 ranks, one run at a time on an otherwise idle machine
(load 2), 3 repeats. Time of the PETSc stage `Second_Run` = one time step with everything compiled and set up.
Minimum of the three repeats; iteration counts of the same runs.

| configuration | iterations | star (s) | column (s) | column / star |
|---|---|---|---|---|
| $K = 100$, $\Delta t = 50$ | 62 → 62 | 16.14 | 11.81 | **0.73** |
| $\Delta t = 5000$ | 29 → 31 | 8.12 | 6.27 | **0.77** |
| $L_z = 250$ m, $\Delta t = 10^5$ | 40 → 44 | 11.23 | 8.94 | **0.80** |
| $\Delta t = 10^5$ | 18 → 20 | 5.30 | 4.35 | **0.82** |
| $N_x = 80$, $\Delta t = 50$ | 11 → 14 | 12.93 | 11.46 | 0.89 |
| $L_z = 12$ km, $\Delta t = 50$ | 13 → 17 | 3.98 | 3.72 | 0.93 |
| $\Delta t = 1000$ | 20 → 29 | 5.82 | 5.96 | 1.02 |
| nominal, $\Delta t = 50$ | 13 → 24 | 4.05 | 5.03 | 1.24 |
| $N_z = 80$, $\Delta t = 50$ | 11 → 35–46 | 13.18 | 26.34 | **2.00** |
| $L_z = 250$ m, $\Delta t = 50$ | 3 → 27–30 | 1.35 | 5.57 | **4.12** |

Where the time goes (nominal, per FGMRES iteration): star 0.31 s, column 0.21 s. The patch solves (`MatSolve`) drop
from 0.125 s to 0.056 s per iteration, a factor 2.2, close to the ratio of the factor sizes (2.5). The rest of the
preconditioner application (residuals of the GMRES smoother, transfers, scatters of the additive Schwarz method:
0.16 s per iteration with the star patch, 0.14 s with the column patch) hardly depends on the patch. Patch factorisation is small in both cases (0.23 s and
0.19 s per time step). So the smaller patch buys about 30 % per iteration, and wins only when it costs fewer than
about 1.4 times the iterations.

**Noise.** The three repeats differ by up to 27 %, in both patches, with two clearly separated values per
configuration (e.g. star nominal 4.05, 4.86, 4.95). I did not find the cause (the machine was idle; CPU frequency
states are my guess, not a finding). Ratios between 0.85 and 1.15 should be read as "no difference".

---

## 6. Conclusions

1. Keep the star patch as the default. It is robust in every parameter, close to an exact solver for field 1,
   and reproducible.
2. The single-column patch is a legitimate option for degree 2 when $\rho_1 = mL_z^2/(\pi^2\Delta x^2)\gtrsim10$
   (large time steps on the nominal geometry: $\Delta t\gtrsim5000$ s) or when the shift is large ($K\gtrsim100$).
   There it is 18 to 27 % faster per time step. Outside that range it is slower, by up to a factor 4.
3. It is not an option for degree 1.
4. If the claim goes into the paper, the supporting argument is Section 2.4 (bubbles plus coarse level), not the
   thin-domain limit. The thin-domain limit is where the column patch is worst.
5. A direction I did not test: the splitting behind the column patch is the hierarchical one in $x$. A smoother
   that adds the vertex functions of the stream function on a *single* vertex column (a Hiptmair-type correction in
   the potential space, $2N_z - 1$ unknowns per column for degree 2) would restore the star patch's constant
   with patches far smaller than 608. Firedrake has `HiptmairPC`; whether it can be used for the mixed
   $(u, u_y, b)$ block on this extruded mesh I have not checked.

---

## 7. How to reproduce

```bash
cd yr2/Boussinesq_implicit
# one run
mpiexec --bind-to none -n 2 python Nonlinear_Boussinesq_Irk_SC.py --nx 20 --nz 20 --length 300000.0 --height 4000.0 \
    --dt 50.0 --tmax 100.0 --shift 1.096622711232e-08 --maxit 100 --dx_test --no_rotation --patch vanka \
    -firedrake_1_snes_atol 1e-09

# the sweeps (each writes ../data_claude_paper_nonlinear_<test>_patch<...>/)
KL="--K_list 0.1 1.0 10.0 100.0 1000.0"
./run_claude_paper.sh -j 12 -T 3600 -m nonlinear -t "dx dz ar dt" -s _patchvanka -x "--patch vanka --snes_atol 1e-9 $KL"
./run_claude_paper.sh -j 12 -T 3600 -m nonlinear -t "dx dz ar dt" -s _patchstar  -x "--snes_atol 1e-9 $KL"
python claude_column_patch/sweeps/compare_patches.py          # tables + figures

# timing, on an idle machine
claude_column_patch/timing/run_timing.sh -n 2 -r 3
python claude_column_patch/timing/parse_timing.py claude_column_patch/timing/runs

# multigrid of field 1 alone
cd claude_column_patch/firedrake_inner_solve && ./run_sweep.sh mytag --K 0.1 1 10 --rhs random bump

# constants of Section 2.5
python claude_column_patch/my_checks/check_two_level_constant.py
```

Finished jobs are skipped on a re-run of the sweeps. Wall time of the sweeps here: 2 h 26 min for all eight.

---

## 8. What is solid and what is not

**How the analysis was produced.** I started six analysis agents in parallel, each to be followed by a referee.
The session hit its usage limit at about 11:15, before any agent had written its report and before any referee had
started. What survives are their scripts and result files. I read the result files, and I re-derived or re-ran the
items marked "mine" below. The items marked "agent" are numerical output of scripts that I have **not** reviewed
line by line.

| Statement | Basis |
|---|---|
| `--patch star` = committed solver | agent (option comparison, residual histories) + mine (220 sweep runs identical to the earlier sweeps) |
| Patch sizes 364 / 608, patch composition | mine (solver view) + agent (enumeration) |
| Parallel patches complete, extra blocks empty | agent's script, run by me with 2 and 4 ranks |
| Field-1 operator = the form of Section 2.1 | agent |
| Divergence matrix independent of geometry | agent |
| Kernel dimensions per patch | agent, two independent implementations agree; formulas derived by me |
| $CG_2(h) = CG_2(2h)\oplus$ bubbles, and the degree-1 failure | mine (argument in 2.4, rank check in my script) + agent (Firedrake prolongation) |
| $C(\rho)$ and its growth with levels | mine (independent script) = agent |
| $\lambda_{min}$ tables of 2.5, mode decoupling, degree-1 law | agent |
| Sweeps, timing, reproducibility table | mine |
| Multigrid alone: $K$ and $\Delta t$ sweeps | agent's runs; $N_x$, $N_z$, $L_z$, smoothing, degree 1: agent's script, run by me |
| Degree 1 in the full solver | mine |

**Not proven.** Section 2.5 gives the constant of one particular splitting, which bounds the additive method. That
the multiplicative V-cycle with a GMRES smoother inherits the same parameter dependence is supported by the
experiments, not derived.

**Not explained.** The growth of the column-patch counts with $N_z$ in the full solver (20 to 39 at
$\Delta t = 50$) does not appear in the multigrid alone (16, 17, 17, 17) and is not predicted by the theory.
Together with the irreproducibility of Section 4.5 it points to the interaction between the inexact field-1 solve
and the outer FGMRES/Newton iteration. I did not investigate further.

**Scope of the experiments.** One test case (Skamarock–Klemp data, amplitude $3\times10^{-4}$, essentially one
vertical mode), $f = 0$, no Eisenstat–Walker, 3 levels, 2 ranks, degree 2. The random right-hand sides of
Section 4.4 are the only evidence for data with many vertical modes. No run with rotation, no linear code.

**Timing.** 20 to 27 % unexplained scatter between repeats (Section 5).

**Housekeeping left to you.** The eight new data folders take 2.2 GB, almost all of it VTK output in
`jobs/*/Nonlinear_slice_imp_ASM*`. I did not delete anything. `claude_column_patch/timing/test_runs/` is an empty
leftover of a test.

---

## 9. Paper section (added 2026-09-30)

In `../Patch-Preconditioner-for-Boussinesq/` (git repository; nothing committed):

| Change | Where |
|---|---|
| New `\subsection{A smaller patch: single-column smoothing}` (label `sec:single-column`, figure `fig:single-column`, table `tab:single-column`) | `ex_article.tex`, inserted at the end of Section 5 immediately before `\subsection{Discussion}` |
| Figure `Image_patch/column_patch_sweeps.png` (four panels: aspect-ratio sweep, mesh sweeps at $\Delta t = 50$ s, shift sweep, collapse against $\rho_1$) | made by `claude_column_patch/paper_figure.py` from the sweep folders |
| The cleveref/ntheorem work-around block of `ex_article.tex` copied into `ex_supplement.tex` | without it the supplement aborts on its first theorem |

Four referee agents (numbers against data, consistency with the paper, mathematics, editorial) reviewed the first draft; their reports are summarised in `claude_column_patch/paper_review_notes/` and led to: the name "single-column patch" throughout (the paper already calls the star patch a "column patch"), an explicit $\delta_{thr} = \Delta t/(2C_M^2)$, the stopping tolerances of the runs, the finest-level $\Delta x$ in $\rho_1$, the two-level/three-level constants (12.3 / 15.7) stated as such, the $N_z$ growth stated as not predicted, the run-to-run variation of the single-column counts stated, "$k = 2$" instead of "$k \ge 2$", and a corrected timing sentence (the column patch is *no faster*, not slower, at $N_x = 80$ and $L_z = 12$ km).

Four one-line edits outside the new subsection keep the paper consistent with it: the caption sentence after Figure `fig:patch-dofs` ("on one star patch highlighted" instead of "on a single column patch"), the sentence "the choice of a vertex-star patch is essential" in Section 3.4 (now: secures robustness for every $k$; for $k = 2$ a single-column patch also meets the kernel condition once the coarse level is counted), the Section 5 roadmap, and the Discussion sentence on the multigrid analysis. **Not changed, for Molin to decide:** Appendix B and Section 3.4 describe the level smoother as six damped Richardson sweeps ($\theta = 0.2$), while `Nonlinear_Boussinesq_Irk_SC.py` (and hence every run of this study, and presumably the paper's own nonlinear runs) uses six GMRES iterations; the new subsection says so explicitly.

The section states only what Sections 2 to 5 of this file support; the two-level constant is quoted as the constant of one splitting, the timing as two-rank minima of three repeats.

**Build hazard.** `latexmkrc` defines a custom dependency that runs `latexmk` on `ex_supplement.tex`, whose `\myexternaldocument{ex_article}` in turn triggers `latexmk` on the article: once the supplement compiles, the two builds call each other without end (several hundred `latexmk` processes within minutes). Build with

```bash
export PATH=$HOME/.TinyTeX/bin/x86_64-linux:$PATH
cd yr2/Patch-Preconditioner-for-Boussinesq
latexmk -norc -pdf -interaction=nonstopmode ex_article.tex      # 12 s, 53 pages
```

(`-norc` ignores `latexmkrc`; the article then only warns that the supplement's `.aux` is absent).

---

## 10. Paper section revised (2026-10-13)

Prompted by the question whether the shift theory needs small $K$: it does ($K \lesssim 1$ is the admissible range of
`cor:dimensional`; the nonlinear outer solve fails at $K \ge 100$ for $\Delta t \ge 50$ s), and the same $K$ sets the
grad–div coefficient $\Delta t/\delta = 2C_M^2/K$ of field 1. So "the single-column patch matches the star for
$K \gtrsim 100$" is true but useless: that regime is outside the theory and the outer iteration is failing there. The
section was rewritten around this; `sec:single-column` in `../Patch-Preconditioner-for-Boussinesq/ex_article.tex`,
pages 40–43 of the 54-page build.

| Change | Detail |
|---|---|
| New paragraph on the shift regime | $K \lesssim 1$ is the grad–div-dominated regime where the kernel condition decides; the test is boundedness as $K \to 0$ |
| New `tab:single-column-block` | block solve alone, random rhs, $K = 0.1 \ldots 10^3$: $k = 2$ star 2,2,2,2,1 / column 20,17,8,5,3; $k = 1$ star 3,3,2,2 / column $>200$,19,13,6 (data: `firedrake_inner_solve/results_sweep_K.csv`, `results_sweep_deg1_K.csv`) |
| Degree-1 full-solver numbers now in the text | 56, 78, 99 against 11, 13, 13 for $N_x = 10, 20, 40$ (`my_checks/degree1_runs/driver.log`) |
| Bounded-in-every-parameter block-solve numbers in the text | $\Delta t$: 17,17,16,9,6; $N_x$: 11–19; $N_z$: 16–17; $L_z$: 16–18 |
| Shift sentences | removed "vanishes once $K \gtrsim 100$ ... grows again for $K = 0.1$"; now: ratio 1.0–3.3 at $K = 1$, 1.4–5.5 at $K = 0.1$, full solver does not yet show the block solve's $K$-independence over the one decade below threshold; agreement at $K \ge 100$ is outside the admissible range |
| $N_z$ growth and run-to-run variation | stated as drawbacks the constant does not predict; $N_z$ growth attributed to the inexact inner solve interacting with Newton–Krylov |
| Timing paragraph | gain 18–23 % at $\Delta t \ge 5000$ s (was "18 to 27 %" including $K = 100$); the $K = 100$ row flagged in the caption as outside the admissible range |
| Closing paragraph | "workable alternative ... drawbacks: larger constant worst in thin weakly stratified domains (up to 10x iterations, 4x time); $N_z$ growth and variability; nothing for $k = 1$; pays off for $\rho_1 \gtrsim 10$; star remains the default" (was "valid, cheaper smoother at large time steps or large shifts") |
| Figure caption (c) | adds the admissible range and that the ten omitted runs are all at $K \ge 100$ |

Unchanged: the figure, the timing table rows, the kernel-count paragraph, the three cross-references outside the
section (lines ~1137, ~2904, ~4003), the Appendix B Richardson/GMRES discrepancy.

Built with `latexmk -norc -pdf -interaction=nonstopmode ex_article.tex` (TinyTeX): no errors, no undefined references,
no overfull boxes in the section. Not committed. The previous version of the file is in the scratchpad of this
session only; `git diff ex_article.tex` in the paper repo shows the change against the committed draft (884524f).

**Open:** one more decade below the threshold ($K = 0.01$, dt sweep, both patches) would settle whether the full
solver inherits the block solve's $K$-independence.
