# Parameter-robustness tests with the analysis-derived shift

Notes accompanying `claude_paper_shift.py`, `Bash/bash_claude_paper.py`,
`plot_claude_paper.py`, `run_claude_paper.sh`, `claude_paper_run_one.sh` and the
generated `test_claude_paper_<mode>_<test>[<suffix>].sh`.  Data and plots live in
`../data_claude_paper_<mode>_<test>[<suffix>]/`.

Written 2026-09-20/21.  The paper is
`Tex_Repo/Patch-Preconditioner-for-Boussinesq/ex_article.tex`; equation and
section labels below refer to it.

---

## 1. What changed relative to `test_dt.sh` / `test_nonlinear_*.sh`

The old sweeps parametrise the shift as $\delta = C_1/\Delta t$ with a fixed
$C_1$.  Section `sec:asymptotic` of the paper (Proposition `prop:asymptotic`,
Corollary `cor:dimensional`) instead gives the admissible shift

$$\delta_{thr} = \frac{\Delta t}{2C_M^2},\qquad
C_M^2 = \max\big(C_h^2,\; m\,C_v^2\big),\qquad
C_h = \frac{L_x}{2\pi},\quad C_v = \frac{L_z}{\pi},\quad
m = 1 + \frac{\Delta t^2 N^2}{4},$$

whose two branches are the two regimes

| regime | condition | $\delta_{thr}$ | scaling |
|---|---|---|---|
| weakly stratified | $4\varepsilon^2 m \le 1$, i.e. $\varepsilon\Delta t N \lesssim 1$ | $2\pi^2\Delta t/L_x^2$ | $\propto \Delta t$ |
| strongly stratified | $4\varepsilon^2 m > 1$, i.e. $\varepsilon\Delta t N \gtrsim 1$ | $\pi^2\Delta t/(2mL_z^2) \to 2\pi^2/(\Delta t N^2 L_z^2)$ | $\propto 1/\Delta t$ |

meeting at $\Delta t = 1/(\varepsilon N) = L_x/(N L_z)$, which is 7500 s on the
nominal $300\,\mathrm{km}\times4\,\mathrm{km}$, $N = 10^{-2}$ geometry.

**All new runs use $\delta = K\,\delta_{thr}(\Delta t, \varepsilon, N)$** with a
dimensionless safety factor $K$; $K = 1$ is the analytical threshold.  Every
member of a sweep therefore gets the shift the analysis prescribes for *its own*
$\Delta t$ and aspect ratio, and the swept quantity is the group
$\varepsilon\Delta t N$ (or, exactly, $4\varepsilon^2 m$).

A useful consequence for the Skamarock–Klemp data (Remark `rem:thin-faster`: the
data contain only pressure modes with vertical structure): with
$\sigma_{01} = \pi^2\Delta t/(mL_z^2)$ from `eq:mode-eigenvalues`, the damping of
the least damped mode present is

$$\rho_{01} = \frac{K}{K + 2/\min(1,\,4\varepsilon^2 m)},$$

a function of $K$ and $4\varepsilon^2 m$ alone.  So at fixed $K$ the iteration
counts of different aspect ratios must **collapse** when plotted against
$4\varepsilon^2 m$ (approximately against $\varepsilon\Delta t N$), and the
paper's operative threshold `eq:C1-threshold-baroclinic` corresponds to
$K_{op} = 1/\min(1, 4\varepsilon^2 m)$ ($\approx 1400$ on the nominal geometry at
small $\Delta t$, $=1$ in the strongly stratified regime).  Both quantities are
implemented in `claude_paper_shift.py` (`rho_01`, `K_operative`) and drawn on
the plots.  What the data support (Sections 4, 5): the collapse in
$4\varepsilon^2 m$, and the ordering of the counts in $K$ at fixed $\Delta t$.
What they do not: $\rho_{01}$ as a quantitative predictor across $\Delta t$ (it
saturates in the strongly stratified regime while the counts fall again, and is
flat in the small-$\Delta t$ corner while the counts grow), and $K_{op}$ as the
location of the failure boundary (Section 4).

The old parametrisation is recovered as $C_1 = \delta\,\Delta t$
(column `C1_equiv` of every `summary_iterations.csv`); the paper's
$C_1 = 10^{-3}$ at $\Delta t = 50$ s is $K = 1823.8$, and the nonlinear
$C_1 = 10^{-5}$ is $K = 18.24$.  Those two values are used for the regression
runs of Section 4.

## 2. Sweeps

Common $\Delta t \in \{1, 10, 50, 100, 500, 10^3, 5\times10^3, 10^4, 5\times10^4, 10^5\}$ s
($\varepsilon\Delta t N$ from $1.3\times10^{-4}$ to 13 on the nominal geometry).

| test | swept | fixed | runs |
|---|---|---|---|
| `dt` (shift scaling) | $K \in \{0.1, 1, 10, 10^2, 10^3, 10^4, 10^5\}$ | nominal geometry, $20\times20$ | 70 |
| `ar` | $L_z \in \{12, 8, 4, 2, 1, 0.5, 0.25\}$ km ($\varepsilon$ from $4\times10^{-2}$ to $8.3\times10^{-4}$) | $K = 1$; linear also $\Delta t = 2\times10^5, 5\times10^5$ | 70 (84) |
| `dx` | $N_x \in \{10, 20, 40, 60, 80\}$ (coarse mesh, $\times4$ after refinement) | $K = 1$ | 50 |
| `dz` | $N_z \in \{10, 20, 40, 60, 80\}$ | $K = 1$ | 50 |

Linear code `LB_time_slice.py` with `--richardson --reordering` as in the old
scripts (rotation is OFF by default there — its `--no_rotation` flag is
`store_false` and would switch rotation *on*); nonlinear code
`Nonlinear_Boussinesq_Irk_SC.py` with `--no_rotation` ($f = 0$, the paper's
setting).  At the time of the paper's nonlinear data (2026-05-20) that script's
`--no_rotation` was also `store_false`, so the old `test_nonlinear_*.sh` runs
were $f = 0$ too; the flag flipped to `store_true` on 2026-06-01.  `maxit = 100`,
two time steps per run, the second one monitored, exactly as before.

## 3. How to run

```bash
cd yr2/Boussinesq_implicit
./run_claude_paper.sh                       # everything: linear+nonlinear, dt ar dx dz, 10 x 2 ranks
./run_claude_paper.sh -m linear -t "dt ar"  # a subset
./run_claude_paper.sh -p -m linear -t dt    # re-plot only
./run_claude_paper.sh -m linear -t ar -s _K100 -x "--K 100"        # side experiment, own folder
./run_claude_paper.sh -m nonlinear -t dx -s _paperK -x "--K 18.24 --dts 50.0 --ew --snes_atol 1e-9"
./run_claude_paper.sh -m linear -t ar -s _paperK -x "--K 1823.8 --rtol 1e-8"
```

Every run executes in its own `jobs/job_<hash>/` (the solver scripts write
`.pvd` and PETSc text files into the cwd, so concurrent runs must not share
one), OpenMPI is given `--bind-to none` (its `-n 2` default is bind-to-core, which
would pin all concurrent runs to cores 0–1), and finished / failed jobs are
skipped on re-runs (`cmd.done` / `cmd.failed`, compared against the command).  The
serial `test_claude_paper_*.sh` files can be run from inside a data directory
like the old `test_*.sh`.

Generator options worth knowing: `--K`, `--K_list`, `--dts`, `--heights`,
`--nxs`, `--nzs`, `--rotation`, `--ew` (Eisenstat–Walker through
`-firedrake_1_snes_ksp_ew`, without editing the script), `--snes_atol`,
`--atol`, `--rtol`, `--shift_override` (a fixed $\delta$ for every run, for
regression against an old fixed-$C_1$ sweep).  `plot_claude_paper.py` reads the solver's `.out` files plus
`jobs/*/status`, so diverged / capped / crashed runs appear as `div` / `>100` /
`crash` in `summary_table.md` and as hatched cells in the iteration maps.

## 4. Results — linear code

Tables are `summary_table.md` in each folder; the numbers below are KSP
iterations $n$ to the solver's stopping criterion.  **Tolerance:** the paper's
linear data were run at `rtol = 1e-8`, the default of `LB_time_slice.py` until
commit 82cd238 (2026-07-08) changed it to `1e-7`; all linear folders were
therefore (re-)run with `--rtol 1e-8` so that they are like-for-like with the
paper (the first pass at `1e-7` is kept in `old_rtol1e-7_run/` sub-folders and
gives exactly one iteration less almost everywhere).  **Counting:** the monitor is
called at iterations $0..n$; the paper's `error_*_Robust_*.png` and
`tab:mode-damping` show the last plotted index $n-1$ (they drop the final
entry), while its `error_shift.png` shows $n+1$ (`len(error)`).  Everything here
is $n$.

**dt / shift scaling** (`data_claude_paper_linear_dt`, 70/70 converged or capped):

| $\Delta t$ | $K_{op}$ | K=0.1 | 1 | 10 | 100 | 1e3 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|---|---|
| 1 | 1406 | 4 | 4 | 4 | 4 | 5 | 9 | 16 |
| 100 | 1125 | 3 | 4 | 4 | 5 | 8 | 17 | 46 |
| 1000 | 54 | 3 | 4 | 5 | 10 | 25 | 68 | >100 |
| 5000 | 2.2 | 4 | 5 | 9 | 18 | 86 | >100 | >100 |
| 10000 | 1 | 4 | 5 | 8 | 16 | 55 | >100 | >100 |
| 100000 | 1 | 3 | 3 | 4 | 5 | 8 | 16 | >100 |

(at `rtol 1e-8`; the `1e-7` pass is one iteration lower almost everywhere)

* $K \le 1$: 3–5 iterations uniformly over five decades of $\Delta t$, both regimes.
* The boundary between bounded and capped counts is **not** a fixed multiple
  of $K_{op}(\Delta t)$: $K_{op}$ falls from 1400 to 1 across the sweep while
  the largest converging $K$ only moves from $10^5$ ($\Delta t \le 500$ s) to
  $10^3$ ($5\times10^3$–$10^4$ s) and back to $10^4$ ($5\times10^4$–$10^5$ s).
  In absolute terms the boundary is close to a constant shift,
  $\delta \approx 10^{-3}\,\mathrm{s\,m^{-2}}$ (converging at
  $1$–$2.5\times10^{-4}$, capped at $1$–$2.5\times10^{-3}$ for every
  $\Delta t \ge 5\times10^3$ s), i.e. the old sweep's $C_1 \approx
  10^{-2}$ boundary at large $\Delta t$.  The worst $\Delta t$ at fixed $K$ is
  the cross-over $5\times10^3$–$10^4$ s, and counts fall again beyond it because
  $\delta \propto 1/\Delta t$ there (`iteration_map_shift.png`).
* The old $C_1$ sweep (`data_linear_dt_deg2`) and the new $K$ sweep lie on one
  $n(\delta)$ curve at every $\Delta t$ once the old shifts are converted to
  $K = \delta/\delta_{thr}$ (e.g. $\Delta t = 10^4$: old $K = 81 \to 15$
  iterations, new $K = 100 \to 15$).
* `error_shift.png` is the counterpart of the paper's `Image_new/error_shift`;
  the curves now peak at the cross-over instead of decaying monotonically, which is
  the reversal of the $\Delta t$-scaling predicted in `sec:asymptotic`.

**ar** (`data_claude_paper_linear_ar`, $K = 1$, 84/84): 3–6 iterations for every
$\varepsilon$ and $\Delta t$ up to $5\times10^5$ s.  The single bump (5–6
iterations) sits at $\Delta t \approx 1/(\varepsilon N)$ for every aspect ratio
and the seven curves collapse to within one iteration in
`error_AR_vs_4eps2m.png` / `error_AR_vs_group.png`, as predicted.  At
$\Delta t = 50$ s all seven aspect ratios converge in the same number of
iterations with overlapping histories: the $\varepsilon$-ordering of the paper's
Fig. `robust-ar` is a property of the old shift ($\rho_{01}$ from 0.86 down to
0.003), and it disappears at $K = 1$ ($\rho_{01} \le 4\times10^{-4}$ for every
$L_z$).

**dx, dz** (`data_claude_paper_linear_dx`, `_dz`, 50/50 each): counts of 3–6
for all $N_x$ at every $\Delta t$, agreeing to within one iteration across the
five meshes (the coarsest $N_x = 10$ occasionally needs one more, and one fewer at
the cross-over); identical for all $N_z$ except $N_z = 10$ at the cross-over (6 vs 5).

**Regression against the paper's figures** (`*_paperK` folders, linear code at
$K = 1823.8$, i.e. the paper's $C_1 = 10^{-3}$ at $\Delta t = 50$):

| | paper's data, $n$ (printed as $n-1$: 18 13 7 5 4 3 3 in `tab:mode-damping`, curves ending at 7 in the dx/dz figures) | new at `rtol 1e-7` | new at `rtol 1e-8` |
|---|---|---|---|
| AR sweep, $\Delta t = 50$, $L_z = 12 \ldots 0.25$ km | 19 14 8 6 5 4 4 | 16 12 7 5 4 4 4 | **19 14 8 6 5 4 4** |
| dx sweep, $\Delta t = 50$, $N_x = 10 \ldots 80$ | 8 8 8 8 8 | 7 7 7 7 7 | **8 8 8 8 8** |
| dz sweep, $\Delta t = 50$, $N_z = 10 \ldots 80$ | 8 8 8 8 8 | 7 7 7 7 7 | **8 8 8 8 8** |

The old and new error histories agree element-wise to better than 1 % over the
common iterations (e.g. $N_x = 20$: 1, 0.289, 0.0375, 1.34e-3, 1.30e-4, 2.2e-5,
2.8e-6 in both), so the solver is unchanged and only the stopping point moves
with the tolerance.

`data_claude_paper_linear_ar_paperK/error_AR_Robust_dt50.0.png` is the
one-to-one counterpart of the paper's `error_AR_Robust_dt50.0`: same counts,
same ordering, same slopes (the only visible difference is the final converged
point, which the paper's plot drops).  The full paper-K AR sweep over $\Delta t$ gives the collapse with a
large dynamic range (3–93 iterations plus six capped runs, peak at
$\varepsilon\Delta t N \approx 0.5$–$1.5$ for every $\varepsilon$;
`error_AR_vs_group.png` there).

## 5. Results — nonlinear code

Same sweeps with `Nonlinear_Boussinesq_Irk_SC.py` ($f = 0$, Newton–Krylov with
the shifted Schur-complement preconditioner on every Newton step, `rtol = 1e-6`,
`ksp_atol = 1e-9`, no Eisenstat–Walker, as the script stands at HEAD).  The
number reported is the **total** number of KSP iterations over all Newton steps
(`snes_ksp_cum[-1]`; the old `error_shift_nonlinear.png` used `len(error)`, which
is the total plus the number of Newton steps).  The AR/dx/dz sweeps carry
`-firedrake_1_snes_atol 1e-9` (see item 5 of Section 6); on the nominal geometry
it never triggers, so the dt sweep, run without it, is consistent with them.

**dt / shift scaling** (`data_claude_paper_nonlinear_dt`, 42 converged, 12 capped, 16 diverged):

| $\Delta t$ | K=0.1 | 1 | 10 | 100 | 1e3 | 1e4 | 1e5 |
|---|---|---|---|---|---|---|---|
| 1 | 4 | 4 | 6 | 10 | 18 | 42 | 131 |
| 100 | 9 | 12 | 21 | 65 | >100 | div | div |
| 1000 | 12 | 20 | 40 | 125 | >100 | div | div |
| 5000 | 16 | 29 | 58 | >100 | >100 | div | div |
| 10000 | 16 | 28 | 60 | >100 | div | div | div |
| 100000 | 12 | 18 | 27 | 55 | >100 | div | div |

Newton steps: 1–2 for $K \le 100$, more when the inner solves stall.  The
totals are larger than the linear counts at the same $K$ (they sum 1–2 Newton
steps, and the two scripts differ in inner solver: GMRES-smoothed multigrid and
FGMRES here, Chebyshev/GMRES with `--richardson` in the linear code, plus
`rtol = 10^{-6}` against $10^{-8}$), and the failure boundary is lower: the
largest converging $K$ is $10^5$ at $\Delta t = 1$ s, $10^3$ at 10 s, $10^2$ at
50–1000 s, 10 at $5\times10^3$–$10^4$ s and $10^2$ at $5\times10^4$–$10^5$ s, so
at $\Delta t = 50$–$100$ s the nonlinear code already fails *below* the paper's
operative threshold $K_{op} \approx 1100$–$1300$.  The shape is the same as in the
linear case: at fixed $K \le 10$ the total peaks at the cross-over
$5\times10^3$–$10^4$ s and falls again beyond it (`snes_ksp_shift.png` =
`error_shift_nonlinear.png`, `iteration_map_shift.png`).  `div` = Newton failed
within 10 iterations (`ConvergenceError`), always with every inner KSP at the
cap.

**ar** (`data_claude_paper_nonlinear_ar`, $K = 1$, 70/70): totals 3–40, 1–3
Newton steps; the peak moves to larger $\Delta t$ as the domain thins and the
seven curves line up against $\varepsilon\Delta t N$ with the maximum at the
cross-over (`snes_ksp_AR_vs_group.png`).  Unlike the linear case the collapse is
in *position* only: the peak height grows as the domain thins (23 total KSP at
$L_z = 12$ km, 29 at 4 km, 32 at 2 km, 37 at 0.5 km, 40 at 0.25 km), so at the
cross-over thinner domains are somewhat slower in the nonlinear code.  Away from
the peak the counts agree to within a few iterations.

**dx, dz** (`data_claude_paper_nonlinear_dx`, `_dz`, $K = 1$, 50/50 each): totals
agree to within ±2 across $N_x$ and across $N_z$ at every $\Delta t$ (4 at
$\Delta t = 1$, 11–13 at 50, 27–31 at the cross-over, 18–19 at $10^5$ s).
`snes_error_dx_dt50.0.png` / `snes_error_dz_dt50.0.png` are the counterparts of
the paper's `Image_nonlinear/snes_error_{dx,dz}_dt50.0`.

**Regression against the paper's nonlinear figures** (folders
`*_paperK*`: the paper's shifts, $f = 0$ and Eisenstat–Walker ON exactly as in
`data_nonlinear_*_deg2_ew`, which are the paper's data — their PNGs are
byte-identical to `Image_nonlinear/*` — plus `snes_atol = 1e-9`, see item 5
below; the $\Delta t = 5000$ set uses `--shift_override 2e-9` because the
paper's $C_1/\Delta t$ is the same for every height whereas $K\delta_{thr}$ is not
for the two thickest domains):

| figure in the paper | paper's data (total KSP) | new |
|---|---|---|
| `snes_error_dx_dt50.0`, $N_x = 10 \ldots 80$ | 25 24 22 22 22 (3 Newton) | 25 22 22 23 22 (3 Newton; $N_x = 10$ identical step by step, 0/3/12/25) |
| `snes_error_dz_dt50.0`, $N_z = 10 \ldots 80$ | 24 24 23 25 25 | 24 22 22 22 22 |
| `snes_error_AR_dt100.0`, $L_z = 12 \ldots 1$ km | 19 20 17 15 11 | 19 20 14 14 10 (and 8, 6 for 0.5, 0.25 km) |
| `snes_error_AR_dt5000.0`, $L_z = 12 \ldots 1$ km | 12 12 10 11 11 | 10 10 10 10 10 (per step 0/2/7/10 vs the paper's 0/2/7/12, 0/2/7/12, 0/2/7/10, 0/2/7/11, 0/2/6/11) |
| `snes_error_shift_dt500.0`, $C_1 = 10^{-3} \ldots 10^{-6}$ | 33 26 21 16 9 9 8 | 28 27 19 12 10 9 7 |

(A first pass of these regression runs was made with rotation ON on the wrong
assumption that the old scripts ran with rotation; it is kept in
`old_rotation_on_run/` sub-folders and gave 24 24 24 22 22 / 29 24 22 22 22 /
19 20 16 14 10 / 9 9 10 10 10 / 32 27 14 13 12 7 7, already within a few
iterations of the paper's data.)

**Eisenstat–Walker variant of the four sweeps** (`*_ew` folders, $f = 0$,
`--ew --snes_atol 1e-9`, $K \le 10^3$ in the dt test): 218 of 220 runs
converged (two diverged at $K = 10^3$, $\Delta t = 5\times10^3$–$10^4$ s).
With EW the solver takes 2–3 Newton steps with cheaper inner solves: totals at
$K \le 1$ are 7–25 (no-EW: 4–29), the $K = 10^2$ column converges at every
$\Delta t$ except the cross-over (77–93 total at $\Delta t = 500$–$1000$ s where
the no-EW runs are capped), and the failure boundary moves up to
$K \approx 10^3$.  The AR sweep again lines up against $\varepsilon\Delta t N$
(peak 18–27); dx/dz totals agree across meshes to ±1 at most $\Delta t$, with two
exceptions worth knowing: at $\Delta t = 1$ s the finer vertical meshes need 25–34
total instead of 9 ($N_z \ge 40$), and at $\Delta t = 50$ s 13 instead of 8 — the
loose first-step tolerance of EW (`rtol0 = 10^{-2}`) occasionally costs an extra
Newton step there.  These folders are the ones to regenerate the paper's
nonlinear figures from if the text keeps the Eisenstat–Walker statement; the
plain folders are the ones to use if the script stays as it is at HEAD.

## 6. Things found on the way that the author should know

1. **The paper compiles again.**  No TeX existed on this machine; a user-space
   TinyTeX is now in `~/.TinyTeX` (PATH added to `~/.bashrc`, `~/.profile`).
   Independently of that, `ex_supplement.tex` lacked the cleveref/ntheorem patch
   that `ex_article.tex` carries, and because `latexmkrc` builds the supplement as
   a dependency (`\myexternaldocument`), the article build was reported as failed
   even though `ex_article.pdf` was produced.  The patch is now in the supplement
   too; `latexmk -pdf ex_article.tex` exits 0 (51 pages, no undefined references).
2. **`--no_rotation` is inverted in `LB_time_slice.py`** (`store_false`): the
   linear runs are $f = 0$ by default, passing the flag switches rotation on.
3. **Eisenstat–Walker is commented out** in `Nonlinear_Boussinesq_Irk_SC.py` at
   HEAD (live at commit 7771540 when the paper's data were made), while the
   paper's appendix says it is used and the paper's nonlinear figures came from
   `data_nonlinear_*_deg2_ew`.  `--ew` in the generator turns it on via the PETSc
   options database.  The main nonlinear sweeps here use the script as it stands
   (no EW); the `*_ew` folders are the same sweeps with EW on.
3b. **`--rtol` default of `LB_time_slice.py` changed from `1e-8` to `1e-7`** (and
   `--maxit` from 150 to 50) in commit 82cd238 of 2026-07-08; the paper's linear
   figures are at `1e-8`.  The nonlinear script's `--no_rotation` flipped from
   `store_false` to `store_true` on 2026-06-01, so the old nonlinear runs were
   $f = 0$ and the current default is rotation ON.
4. **Initial amplitude**: the paper says $\Delta b_0 = 0.01$, both codes use
   $3\times10^{-4}$; the first residual of the old paper runs (0.8947413) equals
   the new one, so the paper's figures were made at $3\times10^{-4}$.
5. **A diagnostic can abort the nonlinear code.**  `ksp_converged_rate` (in the
   script's monitoring options) fits a regression to the residual history; if a
   Newton step's residual is already below `ksp_atol = 1e-9` the KSP stops at
   iteration 0, the fit has one sample and PETSc aborts with
   `LAPACKgetrf INFO = 2`.  It is deterministic and hits thin domains (small
   initial residual) whenever Newton takes one step more than needed, e.g. with
   Eisenstat–Walker.  Work-around used here: `-firedrake_1_snes_atol 1e-9`
   (`--snes_atol 1e-9` in the generator), which stops Newton before such a KSP is
   launched and changes nothing else (a 0-iteration KSP gives a zero update).
   The proper fix is to drop `'ksp_converged_rate': None` from the script.
6. Both scripts print "proportional constant C1 = shift*dt**1.5", which is not
   the paper's $C_1 = \delta\Delta t$; use `C1_equiv` in the CSVs.
7. The paper's `\TODO` about data with a depth-independent pressure component
   (to see the worst-case cross-over at $\varepsilon\Delta t N \sim 1$ directly)
   is not addressed here; it needs a different initial condition in the solver
   scripts.
8. Statements in Section 5 of the paper that the old data already contradict or
   that need re-wording once the new figures go in: (i) "the observed boundary
   between bounded and diverging iteration counts lies near $C_1 = 10^{-2}$":
   in `data_linear_dt_deg2` $C_1 = 1$ converges at $\Delta t = 100$ s in 90
   iterations, at $10^3$ s in 48 and at $10^4$ s in 15 (the paper's own
   `error_shift.png` shows the $C_1 = 1$ curve from $\Delta t = 100$ s on); in the new
   sweep the boundary in $C_1$ terms rises with $\Delta t$, from above
   $C_1 = 2\times10^{-5}$ ($K = 10^5$ still converges) at $\Delta t = 1$ s to
   $C_1 \approx 10^{-2}$, i.e. $\delta \approx 10^{-3}$ s m$^{-2}$, for
   $\Delta t \ge 5\times10^3$ s (Section 4).
   (ii) The captions say FGMRES; the linear runs behind the figures use plain
   GMRES (`--richardson` sets `ksp_type gmres`), FGMRES is the nonlinear code.
   (iii) The dz caption's "more tightly clustered than in the horizontal sweep"
   is not borne out by either the old or the new histories (mean
   $\log_{10}$(max/min) spread at $\Delta t = 50$: dx 0.48, dz 0.47 old; 0.51,
   0.48 new).  (iv) `tab:mode-damping`'s "Iterations" column is the last plotted
   index, one less than the iteration count.
9. Leftovers the author may delete: the `old_*_run/` sub-folders inside the
   `*_paperK*` and linear folders (superseded passes: `rtol 1e-7`, rotation on,
   `atol` variants) and the per-run `jobs/` directories (command, log, status and
   `.out` files of every run; the VTK output was already removed, so all 16
   folders together are about 0.1 GB).
