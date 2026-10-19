# section_referee agent report (workflow wf_d02cd389-216, 2026-10-19)

**Referee report on `\subsection{A smaller patch: single-column smoothing}`** (ex_article.tex:3832-3902), its two figures, and the edits at 1137-1138 and 3925-3927.

No existing file was modified. Everything I wrote is under `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/section_referee/`:
- `check_counts.py` computes every statistic below from `results/all_counts.csv`.
- `paper_fig_single_column_fixed.py` is a copy of the figure script with the cosmetic fixes; its output is in `fig_fixed/single_column_{linear,nonlinear}.png`.
- `latex_build/` holds two pdflatex passes of the current tex.

"Star" below means the 2026-10-16 star re-run. K = δ/δ_thr with δ_thr = Δt/(2C_M²), from Corollary cor:dimensional.

**Verdict.** The section is short and the motivation is close to faithful. But:
- Three statements are false or overstated against the data (1a, 1b, 1c).
- The threshold sentence reads ambiguously and hides that the single-column penalty lies inside the admissible range (2).
- The 30% figure no longer has any support in the paper (1f).
- The captions miss two facts (3).

All other numbers check.

## 1. Numbers and claims

**1a. False: "same counts to within one iteration in all but the thinnest domains" (3857-3861).** In 26 of the 191 linear runs where both patches converge, the counts differ by more than one iteration:
- **Nominal geometry, smallest shift C_1 = 1e-5:** 6, 7, 5 against 4, 4, 3 at Δt = 50, 100, 500 s (`summary_linear_dt.md`).
- **11 thin-domain runs (L_z ≤ 2 km, Δt ≥ 100 s):** including L_z = 2 km at 500 s (6 vs 4) and L_z = 250 m at 100 s (6 vs 4).
- **12 slowly converging runs (Δt ≤ 5 s or C_1 ≥ 0.1):** here the column needs fewer iterations, e.g. 48 vs 62 at N_x = 10, Δt = 1 s; 46 vs 52; 87 vs 90. The only exception is 75 vs 73 on the nominal mesh at Δt = 1 s.

These parts of the sentence are correct:
- "2.75" and "11 against 4 at L_z = 250 m" (Δt = 500 and 1000 s).
- At Δt = 50 s the counts are exact for every mesh (dx and dz: 8 everywhere) and within one iteration for every aspect ratio (19/19, 14/13, 8/8, 6/5, 5/5, 4/5, 4/5).

**1b. False: "the counts coincide where … (Δt ≤ 100 s for C_1 = 1e-5)" (3864-3867).**
- **C_1 = 1e-5, Δt ≤ 10 s:** identical in 43 of 46 runs. The exceptions are 9/8, 28/27, and 30/3 vs 14/2 at N_x = 40, Δt = 10 s.
- **C_1 = 1e-5, Δt = 50-100 s:** identical in only 2 of 32 runs. The column/star median ratio is 1.23 and the maximum 1.80 (18 vs 10 at L_z = 1 km, Δt = 100 s; 25 vs 14 at N_z = 60-80, Δt = 100 s).
- **The paper's own figure contradicts it:** panel (b) of fig:single-column-nonlinear, at Δt = 50 s, shows 24-27 against 22-25.

Over all sweeps the transition is gradual:

| K range | Column/star ratio |
|---|---|
| K ≥ 100 | identical in 53 of 57 runs |
| 10 ≤ K < 100 | median 1.12, max 1.40 |
| 1 ≤ K < 10 | median 1.36, max 2.00 |
| K < 1 | median 2.44, range 1.07-5.88 |

The same overstatement appears in `claude_small_patch_robustness.md:21`. Its own table (lines 172-186) says "identical counts in 58 of the 101 runs" with K ≥ 1.

**1c. Overstated: "its iteration counts are bounded in Δx, Δz, ε, Δt and C_1" (3856-3857).**
- **Vertical refinement, nonlinear, Δt ≥ 500 s:** the column median is 19, 26, 28, 33, 38 for N_z = 10, 20, 40, 60, 80; the star stays at 10. The extreme is 47 vs 8 at N_z = 80, Δt = 5e4 s. No figure shows this.
- **Aspect ratio (ε):**
  - Linear, Δt = 500-1000 s: the column goes from 5 to 11 between L_z = 4 km and 250 m, while the star goes from 5 to 4.
  - Nonlinear, Δt = 500 s: the column goes from 14 at L_z = 12 km to 37 at 1-2 km; the star stays at 9-12.
- **Δt and C_1:** neither patch is "bounded" in Δt or C_1. Counts range from 4 to 90 across Δt, and large C_1 hits the iteration cap. The paper's own claim (3304-3306) is "robust in Δx, Δz, ε, and (within the analytically predicted range of) C_1".
- **What is true:** at Δt = 50 s, the time step of the mesh sweeps, the column counts are as flat as the star's. Linear counts are identical. Nonlinear: N_x 24-26 vs 22-25; N_z 24-27 vs 22-24; L_z 14-27 vs 10-24.

**1d. Correct with a qualifier: "two to three times" (3867-3868).** For C_1 = 1e-5, Δt ≥ 500 s (80 runs, none equal): median 2.60 (per sweep 2.30 dx, 2.44 ar, 2.89 dz), range 1.40-5.88. The text should give the median and the maximum.

**1e. Correct.**
- **"Same set of non-converging runs":** the nonlinear problem has the same 13 failures for both patches (ar: L_z = 12, 8 km at 1 s; dt sweep: 11).
- **"Converges in every run in which the star converges":** in the linear problem the column also converges at N_z = 40, Δt = 1 s (96 iterations) where the star hits the cap; failures are 15 vs 16.
- **"364 against 608":** 18N_z+4 and 30N_z+8 follow from the spaces:
  - u: 2 facets × 2N_z plus the bubble's 2N_z.
  - w and b: 4N_z+2 each per column.
  - u_y: 4N_z per column.

  They are consistent with Firedrake's ASMVankaPC (closure of star; `asm.py:191-197`), as used at `LB_time_slice.py:230-236` and `Nonlinear_Boussinesq_Irk_SC.py:210-216`.

**1f. No support left in the paper: "about 30% … factor of about 1.4" (3870-3872).**
- 1.4 = 1/0.7, so the two numbers are consistent with each other.
- The 30% comes from `claude_column_patch_study.md:41` and `:498-503`. That was the nonlinear solver with a GMRES(6) smoother at (20,20) on 2 ranks: 0.31 vs 0.21 s per FGMRES iteration, with patch solves 2.2× cheaper.
- The table that showed this was deleted today, and the linear solver (Chebyshev(6) smoother) was never timed. The text should state where the figure comes from.
- The text should also say where the counts are within the 1.4 factor:
  - **Linear:** 177 of 191 runs. It fails in the 11 thin-domain runs and the 3 runs at C_1 = 1e-5, Δt = 50-500 s.
  - **Nonlinear:** 93 of the 101 runs with K ≥ 1, but only 7 of the 106 with K < 1.

**1g. "With the same … tolerances" (3853-3854).** This is true between the two patches. But the nonlinear re-runs add `-firedrake_1_snes_atol 1e-9` (`data_nonlinear_*_ew_*/commands.txt`; `claude_paper_tests_README.md:317-325`). So the star curves in the nonlinear figure are not the counts behind Figs. snes-* (e.g. Δt = 50 s, N_z = 20-80: 22 vs the paper's 23-25). Optional caption clause in R6.

## 2. Direction of the threshold

The admissible condition is an upper bound on the shift: δ ≲ Δt/(2C_M²). It appears at eq:delta-admissible (2679-2682) and Corollary cor:dimensional (2765-2770). For these runs it is eq:C1-threshold (3216-3221): C_1 ≲ 2π²Δt²/L_x² ≈ 2.2e-10 Δt².

So K ≥ 1 means the shift is at or above the bound, i.e. outside the admissible range; K < 1 is inside it.
- For C_1 = 1e-5 the crossover is at Δt ≈ 214 s, for every member of the sweeps.
- All L_z are weakly stratified at Δt ≤ 100 s.
- With L_z = 12 km at large Δt, the strongly stratified branch requires C_1 ≲ 1.4e-3, which still holds.
- K = 18 at 50 s, 4.6 at 100 s, 0.18 at 500 s.

The numbers "(Δt ≤ 100 s for C_1 = 1e-5)" and "below it (Δt ≥ 500 s)" are therefore right, and the usage matches line 3620 ("operates above the worst-case threshold"). But "at or above the admissible threshold" invites the opposite reading. It also hides that the 2-3× penalty is inside the admissible range, the regime the analysis recommends. The sentence should say "outside" and "inside" explicitly (R4).

## 3. Captions against figures and `paper_fig_single_column.py`

**Correct:**
- Labelled values (`:46`, `:72`): linear L_z 12, 2, 0.25 km of 7 and C_1 = 1, 1e-3, 1e-5 of 6; nonlinear L_z 12, 4, 1 km of 5 and C_1 = 1e-3, 5e-5, 1e-6 of 7.
- Colour ramp: viridis 0→0.82, dark = 12 km.
- Line and marker styles: star solid with filled circles, column dashed with open squares; in panel (b) circles for N_x and squares for N_z.
- Linear (b), "four curves coincide at 8": all 20 counts are 8.

**Linear caption: "Runs that reached the iteration cap are omitted".** The script drops every status ≠ 'ok' (`:49`, `:75`). That also drops the DIVERGED_BREAKDOWN run at C_1 = 1, Δt = 1 s. Panel (b) is not filtered (`:63`), but every Δt = 50 s run converged. In total 10 runs are omitted per patch, the same set for both.

**Linear (c):** the colour ramp is not described.

**Nonlinear caption:** panel (b) lacks "C_1 = 10^{-5}", and there is no omission note although 13 runs per patch are dropped (2 in (a), 11 in (c)).

**Cosmetic problems** (visible on p. 42 of the compiled PDF):
- In panel (b) of both figures the x-tick labels 10 and 20 collide and print as "1020".
- The panel (b) legend (ncol=2) runs past the right spine towards panel (c).
- Panel (b) has fractional y ticks.
- The legend at 6.5 pt prints at about 5 pt (scale 370/474).

Fixes, checked in `fig_fixed/`:
- Line 36: `fontsize=7.5` (or 8).
- Line 65: `ax.set_xscale('log'); ax.set_xticks([10, 20, 40, 60, 80]); ax.set_xticklabels(['10', '20', '40', '60', '80']); ax.xaxis.set_minor_locator(plt.NullLocator())`, then add `ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))`.
- Line 67: `ax.legend(loc='lower right', ncol=1, **LEG)`.

## 4. Motivation and length

Lines 3836-3841 give Molin's reason briefly and as motivation ("suggested", "expected"): a thin domain makes the divergence mostly a vertical derivative, so a single column should hold what the patch needs, hence the Vanka patch with construct_dim 1 on mg_levels. Three problems:
- **Misstated condition:** "the divergence-free fields that a robust patch must contain" does not match the condition at 1141-1144, which asks for a basis of the kernel whose members are each supported on one patch.
- **Results contradict the motivation silently:** the single column does worst in the thinnest domains. Linear: within one iteration for L_z ≥ 4 km, 2.75× at 250 m. Nonlinear at Δt = 500 s: 37 at 1-2 km against 14 at 12 km.
- **The k = 1 sentence also contradicts it**, since for k = 1 a column holds no divergence-free field at all.

Do not add Molin's literal "the kernel will be smaller". The discrete kernel does not depend on ε (`claude_column_patch_study.md:16-18`, and the committed HEAD text of this section). Two paragraphs is the right length; R1 plus R4 add about 14 lines.

## 5. LaTeX and consistency

- **Compiles:** no undefined references or citations, no overfull boxes in 3832-3902, and both figures land on p. 42.
- **Labels:** all resolve (sec:inner-solvers 1076, app:petsc-options 4262, cor:dimensional 2766, sec:robustness 3027, sec:nonlinear 3320).
- **`\opt`:** defined as a URL command (line 77), and the colon style matches 4265-4266.
- **Other edited sentences:**
  - Lines 1137-1138 are consistent with the section.
  - Lines 2903-2904 are consistent.
  - The Discussion sentence (3925-3927) is grammatical. It writes "grad-div" where the rest of the paper uses "grad--div"; this predates today's edit.
- **Terminology:** the paper calls the star patch the "column patch" or "column-patch smoother" at 1123, 1284, 1291, 1298, 3715, 3725-3726, 3922, 4089, 4116, 4261 and 4264. "Single-column patch" next to "column patch" can confuse; R3 ties them together.
- **Coverage:** the single-column patch is not mentioned in the abstract, introduction or conclusions.

## 6. Replacement text

**R1, lines 3836-3841**

Old:
```
A smaller patch is
suggested by the geometry: as the aspect ratio decreases the domain
becomes thinner and the divergence is dominated by its vertical
derivative, so the divergence-free fields that a robust patch must
contain (Section~\ref{sec:inner-solvers}) are expected to be captured
within a single column of cells.
```
New:
```
A smaller patch is
suggested by the geometry.  A robust patch smoother needs the
divergence-free fields, the kernel of the grad--div term, to have a
basis of functions each supported in one patch
(Section~\ref{sec:inner-solvers}); as the aspect ratio decreases the
domain becomes thinner and the divergence is dominated by its vertical
derivative, so this kernel is expected to be captured within single
columns of cells.
```

**R2, lines 3846-3847**

Old:
```
\opt{pc_vanka_construct_dim:1} in place of \opt{firedrake.ASMStarPC}
(Appendix~\ref{app:petsc-options}).
```
New:
```
\opt{pc_vanka_construct_dim:1} in place of \opt{firedrake.ASMStarPC}
and \opt{pc_star_construct_dim:0} (Appendix~\ref{app:petsc-options}).
```

**R3 (optional), lines 3835-3836**

Old:
```
The star patch of Section~\ref{sec:inner-solvers} spans the two
columns of cells around a vertex of the base mesh.
```
New:
```
The star patch of Section~\ref{sec:inner-solvers}, the column patch
of the preceding sections, spans the two columns of cells around a
vertex of the base mesh.
```

**R4, replace all of lines 3851-3874.** Every number below was verified. For brevity the clauses "so that the thin domains … case" and "where the grad--div term dominates" can be dropped. Keep the N_z clause if "parameter-robust" stays.
```
Figures~\ref{fig:single-column-linear} and~\ref{fig:single-column-nonlinear}
show the four sweeps of Sections~\ref{sec:robustness}
and~\ref{sec:nonlinear}, with the same shifts $\delta = C_1/\Delta t$,
meshes, aspect ratios and tolerances, for both patches.  The
single-column patch converges in every run in which the star patch
converges, and at $\Delta t = 50\,\mathrm{s}$, the time step of the mesh
sweeps, its counts depend as little on $\Delta x$, $\Delta z$ and
$\varepsilon$ as those of the star patch.  For the linear problem the
two patches give the same counts to within one iteration in $165$ of
the $191$ runs in which both converge, among them every run of
Figures~\ref{fig:robust-dx}--\ref{fig:robust-ar}; the single column
needs clearly more only for the smallest shift, $C_1 = 10^{-5}$ (at
most $7$ against $4$), and in thin domains at large time steps
($L_z \le 2\,\mathrm{km}$, $\Delta t \ge 100\,\mathrm{s}$), up to $2.75$
times the iterations ($11$ against $4$ at $L_z = 250\,\mathrm{m}$,
$\Delta t = 500$--$1000\,\mathrm{s}$), so that the thin domains that
suggested it are its least favourable case.  For the nonlinear problem
the two patches fail in the same runs.  Their counts agree where the
shift lies far above the worst-case admissible bound
$\delta \lesssim \Delta t/(2C_M^2)$ of Corollary~\ref{cor:dimensional},
that is, well outside the admissible range ($\Delta t \le 10\,\mathrm{s}$
for $C_1 = 10^{-5}$), and differ by at most a factor $1.8$ at
$\Delta t = 50$--$100\,\mathrm{s}$, still above the bound.  Inside the
admissible range ($\Delta t \ge 500\,\mathrm{s}$), where the grad--div
term dominates, the single column needs two to three times the
iterations (median $2.6$, at most $5.9$), and its counts grow under
vertical refinement (median $19$ at $N_z = 10$, $38$ at $N_z = 80$,
against $10$ for the star).  The smaller patch is therefore a working,
parameter-robust alternative to the star patch, at about $30\%$ lower
cost per outer iteration (measured for the nonlinear solver at
$(N_x, N_z) = (20, 20)$), so that it is also faster wherever the two
counts are within a factor of about $1.4$: in the linear problem except
for $C_1 = 10^{-5}$ at $\Delta t = 50$--$500\,\mathrm{s}$ and in thin
domains, and in the nonlinear problem outside the admissible range.
All runs use $k = 2$; for $k = 1$ a single column contains no
divergence-free fields and the star patch is required.
```

**R5, lines 3886-3888**

Old:
```
(c) Shift sweep on the nominal
    geometry; three of the six values of $C_1$ are labelled.  Runs
    that reached the iteration cap are omitted.}
```
New:
```
(c) Shift sweep on the nominal
    geometry; three of the six values of $C_1$ are labelled, the colour
    runs from $C_1 = 1$ (dark) to $10^{-5}$ (light).  Runs that did not
    converge within the iteration cap, the same ten for both patches,
    are omitted.}
```

**R6, lines 3898-3900** (the last sentence is optional; see 1g)

Old:
```
(b) Mesh
    resolution at $\Delta t = 50\,\mathrm{s}$.  (c) Shift sweep on the
    nominal geometry; three of the seven values of $C_1$ are labelled.}
```
New:
```
(b) Mesh
    resolution at $\Delta t = 50\,\mathrm{s}$, $C_1 = 10^{-5}$.  (c) Shift
    sweep on the nominal geometry; three of the seven values of $C_1$
    are labelled.  Runs that did not converge, the same thirteen for
    both patches, are omitted.  The runs add a Newton absolute
    tolerance of $10^{-9}$ to the settings of
    Section~\ref{sec:nonlinear}.}
```

**R7 (optional), last sentence of the paragraph, after R4**

Old:
```
All runs use $k = 2$; for
$k = 1$ a single column contains no divergence-free fields and the
star patch is required.
```
New:
```
All runs use $k = 2$, for which a single column contains the
divergence-free fields built on the quadratic bubbles of its cells; for
$k = 1$ it contains none, and the star patch is required.
```

**R8 (optional), line 3926:** `grad-div ` → `grad--div `.

## 7. Outside this section, for the orchestrator

- **Table 5 (tab:monolithic):**
  - Line 3827: the monolithic patch DOFs at (40,80) are printed as 3408; 38N_z+8 gives 3048.
  - Caption, lines 3813-3814: says "ε ∼ 10^{-3}", but the runs use L_z/L_x = 4/300 = 1.3e-2.
- **Appendix app:petsc-options (4271-4316):**
  - It shows the multigrid levels using a Richardson smoother. The runs actually use Chebyshev(6) for the linear problem (`LB_time_slice.py:284`, under `--richardson`) and GMRES(6) for the nonlinear one (`Nonlinear_Boussinesq_Irk_SC.py:250`).
  - The linear outer solver is GMRES (`LB_time_slice.py:330`), while the paper says FGMRES.
- **Run-to-run variability:** the nominal configuration appears in all four sweeps with identical commands.
  - Star counts are identical across the four at every Δt, in both codes.
  - Column counts vary: linear 26/27 at Δt = 5 s; nonlinear 26/28 at Δt = 500 and 5000 s, and 22/24/26 at 5e4 s.
  - Single-run ratios such as 5.9 therefore carry roughly ±10%.
- **Auxiliary-preconditioner defect:** the nonlinear counts of both patches include the T^{-1}A_δ^{-1} defect (`claude_small_patch_robustness.md:28-40`). With the fix, the column/star median in the K = 1 sweep fell from 1.67 to 1.42, so the nonlinear penalty of the single column may be overstated.
