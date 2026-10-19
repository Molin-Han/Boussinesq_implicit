# figure_referee agent report (workflow wf_d02cd389-216, 2026-10-19)

Figure referee report for the ten restyled figures (Figs. 4–12, `sec:robustness` and `sec:nonlinear`). All tex paths below are in `/home/mh2520/molin/phd/yr2/Patch-Preconditioner-for-Boussinesq/`. No paper or solver file was modified.

## 1. The data are the same in all ten figures

- **Restyled files are current.** Each `Image_restyled/*.png` is pixel-identical to a fresh run of `claude_small_patch/restyle_paper_figures.py --set paper` (written to my folder, `regen/`).
- **Originals reproduced from the same data.** I re-drew every figure in the old matplotlib style from the data the restyle script loads, then compared with `Image_new/` and `Image_nonlinear/`.
  - Six figures (Figs. 4, 5, 8, 9, 10 left, 11): 0 differing pixels inside the axes outside the legend box.
  - Four tight-cropped figures (Figs. 6, 7, 10 right, 12): same marker count per colour, centroid shift at most 1.41 px. That is crop rounding: one iteration is at least 4.8 px.
  - Evidence: `replica_check2.out`, `marker_check.out`.
- **End points match.** Figs. 4–5: every run 8 iterations (curves end at index 7). Fig. 6: curves end at 18, 13, 7, 5, 4, 3, 3. Figs. 8/9: totals 25, 24, 22, 22, 22 and 24, 24, 23, 25, 25. Fig. 10: left 19, 20, 17, 15, 11; right 7, 7, 7, 7, 6. Fig. 11: 33, 26, 21, 16, 9, 9, 8. Full listing in `dump_figure_data.out`.
- **Axis ranges are the same**, except that Figs. 7 and 12 now start at y = 0 (harmless). Plot titles were removed and x-labels added; the captions still fit.
- **Legend renames:**
  - nx/nz → N_x/N_z
  - dx=30000.0 / dz=400.0 → N_x/N_z
  - K → C_1 (this now matches the text)
  - AR → L_z in Figs. 6 and 10

  The text never used the old names. It mentions no colours or line styles in these two sections.

## 2. Problems introduced by the restyle (fix before pushing)

**2a. The legend now hides data in three figures.** The originals hid no markers (`legend_overlap.out`). Candidate fixes are in `cand/`, each checked with the same overlap test:

- **Fig. 6** (`error_AR_Robust_dt50.0`): hides L_z = 12 km at iterations 8–12 and L_z = 8 km at iteration 8.
  - Candidate `cand/error_AR_Robust_dt50.0__shortlabels.png`: legend shows L_z only.
  - In the restyle script: `restyle_paper_figures.py:95` `short=(size == ps.HALF)` → `short=(size == ps.HALF or test == 'ar')`.
  - If adopted, apply edit A1 (caption) below.
- **Fig. 7** (`error_shift`): hides the C_1 = 1 point at Δt = 10³ (49).
  - Candidate `cand/error_shift__ncol3_ylim120.png`: 3-column legend, y-axis to 120. In the script, add after line 126 `ax.legend(loc='upper right', ncol=3, columnspacing=1.0)` and change line 127 to `ax.set_ylim(0, 120)`.
  - Same fix with true counts: `cand/error_shift__corrected_ncol3_ylim120.png` (see section 4).
- **Fig. 11** (`snes_error_shift_dt500.0`): hides the C_1 = 10⁻³ Newton-step kink at (27, 1.9e-8).
  - Candidate `cand/snes_error_shift_dt500.0__nolegend.png`: no legend; Fig. 12 beside it uses identical colours for the same seven C_1.
  - If adopted, end the new caption of edit B10 with `...10^{-6}\}$, coloured as in Figure~\ref{fig:snes-shift-sweep}.}`.
  - Alternative: `fix_trial/snes_error_shift_dt500.0__outside_right.png`, which makes the panel wider.

**2b. Fig. 10 legend.** It now shows L_z, not ε. Edit B8 adds the L_z ↔ ε mapping to the caption.

**2c. Images not in git.** `Image_restyled/` is untracked. The tex now needs 12 files from it, so a push without them will not build. Run:

```
git add ex_article.tex Image_restyled/{error_dx_Robust_dt50.0,error_dz_Robust_dt50.0,error_AR_Robust_dt50.0,error_shift,snes_error_dx_dt50.0,snes_error_dz_dt50.0,snes_error_AR_dt100.0,snes_error_AR_dt5000.0,snes_error_shift_dt500.0,error_shift_nonlinear,single_column_linear,single_column_nonlinear}.png
```

Add `Image_restyled/corrected/error_shift_nonlinear.png` too if edit C1 is adopted, plus any adopted candidates.

**2d. Build log** (`ex_article.log`, 06:13, 53 pages):
- All ten images are found, at the expected lines (3046, 3063, 3187, 3289, 3348, 3356, 3375, 3379, 3409, 3417) on pp. 30, 30, 32, 34, 35, 35, 35, 35, 36, 36.
- No missing files, no undefined references or citations, no rerun request.
- No warning comes from the image swap itself, with one possible exception: two `Underfull \vbox (badness 10000)` while shipping page 35 (Figs. 8–10). This may come from the taller panels; every restyled figure is taller at the same width, e.g. 142.6 pt vs 133.3 pt for the nonlinear panels. Page 35 looks normal.
- Other warnings are in unchanged text: overfull boxes at tex 3017–3025 (55 pt), 3753–3766 (71 pt) and 3817–3830 (26 pt), and the hyperref bookmark and cmtt font warnings.
- The log predates the last tex edit (06:45): it places `single_column_*` at lines 3874/3890, the file now has 3878/3894. Rebuild before committing, with `latexmk -norc -pdf -interaction=nonstopmode ex_article.tex`.

## 3. Sentences the plotted data contradict

These were already wrong with the original figures, which show the same data; they count as "counts read off the figure". Exact strings, trailing spaces included, are in `edits.json`. Each old string was checked to occur exactly once and all edits were dry-run on a copy (`ex_article.figure_referee.tex`, `edits.diff`); braces stay balanced. I did not compile the copy, because the timing runs are using the cores.

**B1 — Fig. 5 text, tex 3057–3059.** Every run takes 8 iterations in both sweeps. The spread across curves is not smaller for N_z: 2.52× vs 2.32× at the last plotted point, and 7.9× vs 3.5× at iteration 3.

OLD:
```
The iteration counts are even more tightly clustered than in the 
horizontal sweep, reflecting the fact that the column-patch 
smoother resolves the vertical coupling exactly within each patch.
```
NEW:
```
As in the horizontal sweep, every run takes the same number of
iterations, reflecting the fact that the column-patch
smoother resolves the vertical coupling exactly within each patch.
```

**B2 — Fig. 5 caption, tex 3067–3070.**

OLD:
```
    held at their nominal values. The curves are more tightly
    clustered than in Figure~\ref{fig:robust-dx}, reflecting the
    fact that the column-patch smoother resolves the vertical
    coupling exactly within each patch.}
```
NEW:
```
    held at their nominal values. As in
    Figure~\ref{fig:robust-dx}, every curve ends at the same
    iteration, reflecting the fact that the column-patch smoother
    resolves the vertical coupling exactly within each patch.}
```

**B3 — Fig. 7 text, tex 3226–3229.** C_1 = 1 converges at Δt = 100 s (90 iterations, plotted as 91). C_1 = 10⁻² reaches the cap at Δt = 1 s.

OLD:
```
iteration count diverges (the $C_1 = 1$ curve never converges to the
tolerance within the iteration cap for $\Delta t \le 10^2$), while for
$C_1 \le 10^{-2}$ it remains bounded and decreases as $\Delta t$
grows;
```
NEW:
```
iteration count diverges (the $C_1 = 1$ curve reaches the iteration cap
for every $\Delta t \le 50\,\mathrm{s}$, and the $C_1 = 10^{-1}$ curve for
every $\Delta t \le 10\,\mathrm{s}$), while for $C_1 \le 10^{-3}$, and for
$C_1 = 10^{-2}$ when $\Delta t \ge 5\,\mathrm{s}$, it remains bounded and
decreases as $\Delta t$
grows;
```

**B4 — Fig. 7 text, tex 3258.** At Δt = 5 s, C_1 = 10⁻⁵ needs 7 (plotted), against 5 at Δt ≥ 50 s, so "large" does not hold.

OLD: `counts are large for every $C_1$ and ordered monotonically in it.`

NEW:
```
counts exceed their values at $\Delta t \ge 50\,\mathrm{s}$ for every
$C_1$ and increase monotonically with it.
```

**B5 — Fig. 7 caption, tex 3298–3300.** The C_1 = 1 and 10⁻¹ curves do not level off, and the omitted capped runs are not mentioned.

OLD:
```
    near $10^{-2}$, for $\Delta t \gg 2/N$.  The curves separate
    accordingly into those that respect it and those that do not, and
    level off beyond $\Delta t \sim 2/N$.}
```
NEW:
```
    near $10^{-2}$, for $\Delta t \gg 2/N$.  The curves separate
    accordingly into those that respect it, which level off beyond
    $\Delta t \sim 2/N$, and those that do not.  Runs that reach the
    iteration cap of $100$ are omitted.}
```

**B6 — tex 3330.** Figs. 8–10 used C_1 = 10⁻⁵ (`Bash/bash_nonlinear.py:14,28,35`). Their histories are identical to the C_1 = 10⁻⁵ runs of the Δt sweep at Δt = 50, 500 and 5000 s. The text at 2996–2997 says the nominal C_1 = 10⁻³ applies "throughout".

OLD: `preconditioned by $A_\delta^{-1}$ as before. The convergence`

NEW:
```
preconditioned by $A_\delta^{-1}$ as before, with
$\delta = C_1/\Delta t$ and $C_1 = 10^{-5}$ in the resolution and
aspect-ratio sweeps of Figures~\ref{fig:snes-dx}--\ref{fig:snes-ar}. The convergence
```

**B7 — Fig. 10 text, tex 3367–3369.** At Δt = 100 s the totals are 19, 20, 17, 15, 11, which is not "tightly clustered".

OLD:
```
$\Delta t = 5000\,\mathrm{s}$. In both cases the curves are
tightly clustered, confirming that the preconditioner retains
its aspect-ratio robustness in the nonlinear setting.
```
NEW:
```
$\Delta t = 5000\,\mathrm{s}$. At $\Delta t = 5000\,\mathrm{s}$ the
curves are tightly clustered; at $\Delta t = 100\,\mathrm{s}$ the
cumulative count falls from $19$--$20$ at $L_z = 12$ and $8\,\mathrm{km}$
to $11$ at $L_z = 1\,\mathrm{km}$, the same improvement in thinner
domains as in the linear sweep (Remark~\ref{rem:thin-faster}).  In
neither case does the count grow as $\varepsilon$ decreases, confirming
that the preconditioner retains its aspect-ratio robustness in the
nonlinear setting.
```

**B8 — Fig. 10 caption, tex 3384–3386.** The right panel drops the converged Newton step and the left panel keeps it; the true right-panel totals are 12, 12, 10, 11, 11.

OLD:
```
    (left) and $\Delta t = 5000\,\mathrm{s}$ (right) for five
    aspect ratios. The curves are tightly clustered at both time
    steps.}
```
NEW:
```
    (left) and $\Delta t = 5000\,\mathrm{s}$ (right) for five
    aspect ratios, $L_z \in \{12, 8, 4, 2, 1\}\,\mathrm{km}$
    ($\varepsilon = L_z/L_x$ from $4\times10^{-2}$ to $3.3\times10^{-3}$).
    Thinner domains converge in fewer iterations at
    $\Delta t = 100\,\mathrm{s}$ ($11$ against $19$--$20$ cumulative
    iterations); at $\Delta t = 5000\,\mathrm{s}$ the curves are tightly
    clustered.  In the right panel the final Newton step, which takes
    the error below $10^{-12}$, is not plotted.}
```
Option: use `cand/snes_error_AR_dt5000.0__keepconverged.png` for the right panel (keeps the converged point, no legend overlap; script line 212, drop the trailing `True`). Then replace the last two sentences with `at $\Delta t = 5000\,\mathrm{s}$ the curves are tightly clustered ($10$--$12$ cumulative iterations).}`.

**B9 — Figs. 11–12 text, tex 3397–3398.**
- C_1 ∈ [10⁻⁴, 10⁻³] gives the highest counts, not uniformly low ones:
  - C_1 = 10⁻⁴ needs 74–90 at Δt = 5–50 s.
  - C_1 ≥ 5×10⁻⁴ exceeds 100 at Δt ≤ 100 s.
- At Δt ≤ 500 s the count falls monotonically as C_1 decreases, and C_1 = 10⁻⁶ is lowest.

OLD:
```
of Section~\ref{sec:robustness}: a well-chosen $C_1$ in the range
$[10^{-4}, 10^{-3}]$ gives uniformly low iteration counts, and
```
NEW:
```
of Section~\ref{sec:robustness}: for $\Delta t \le 500\,\mathrm{s}$ the
iteration count decreases monotonically as $C_1$ decreases, so that the
smallest value tested, $C_1 = 10^{-6}$, gives the lowest counts, while
for $C_1 \ge 5\times10^{-4}$ and $\Delta t \le 100\,\mathrm{s}$ the count
exceeds $100$; for $\Delta t \ge 10^{4}\,\mathrm{s}$ all seven values give
nearly the same count; and
```

**B10 — Fig. 11 caption, tex 3410–3412.** The interval is reversed, and the values are seven discrete ones.

OLD:
```
        \caption{Nonlinear shift-parameter sensitivity: convergence
        history at $\Delta t = 500\,\mathrm{s}$ for
        $C_1 \in [10^{-3}, 10^{-6}]$.}
```
NEW:
```
        \caption{Nonlinear shift-parameter sensitivity: convergence
        history at $\Delta t = 500\,\mathrm{s}$ for the seven values
        $C_1 \in \{10^{-3}, 5\times10^{-4}, 10^{-4}, 5\times10^{-5},
        10^{-5}, 5\times10^{-6}, 10^{-6}\}$.}
```

**B11 — Fig. 12 caption, tex 3418–3420.** This text assumes edit C1 is applied.

OLD:
```
        \caption{Nonlinear shift-parameter sweep: total KSP
        iteration count vs $\Delta t$ for
        $C_1 \in [10^{-3}, 10^{-6}]$.}
```
NEW:
```
        \caption{Nonlinear shift-parameter sweep: total KSP
        iteration count over the Newton steps vs $\Delta t$ for the
        seven values of $C_1$ of Figure~\ref{fig:snes-shift-hist}.
        Runs in which a Newton step reached the cap of $100$ KSP
        iterations, or which did not complete, are omitted.}
```

## 4. Iteration counts are plotted inconsistently

The restyle keeps the old plotting rules on purpose. As a result:
- **Figs. 4–6** end each curve one index before the true count n (the converged point, error below 1e-13, is dropped). The 0.5 km curve therefore stops at 1.7e-5.
- **Fig. 7** plots n+1.
- **Fig. 12 (drop-in)** plots total KSP iterations plus the number of Newton steps. At Δt = 500 s it shows 36, 29, 24, 19, 11, 11, 10, while Fig. 11 next to it ends at 33, 26, 21, 16, 9, 9, 8. It also omits C_1 = 5×10⁻⁴ at Δt = 100 s, which converged in 129 iterations.
- **`sec:single-column` figures** use true counts.

So the nominal run (Δt = 50 s, C_1 = 10⁻³) appears as 7 (Figs. 4–6, `tab:mode-damping`), 9 (Fig. 7) and 8 (`fig:single-column-linear` caption, about tex 3885).

Recommendation:
- **C1, tex 3417:** `{Image_restyled/error_shift_nonlinear}` → `{Image_restyled/corrected/error_shift_nonlinear}`. That figure shows true totals, has no legend overlap, and is current.
  - If C1 is not adopted, use this B11 text instead: `\caption{Nonlinear shift-parameter sweep: total KSP iterations plus the number of Newton steps vs $\Delta t$ for the seven values of $C_1$ of Figure~\ref{fig:snes-shift-hist}; values of $100$ or more are not plotted.}`
- **Fig. 7:** use `cand/error_shift__corrected_ncol3_ylim120.png` (true n; 8 at the nominal point). If the n+1 version is kept, append ` The value plotted is the iteration count plus one.` to the caption.
- **C2, tex 3050–3051:**

  OLD:
  ```
      fixed $\Delta t = 50\,\mathrm{s}$, with remaining parameters
      held at their nominal values.}
  ```
  NEW:
  ```
      fixed $\Delta t = 50\,\mathrm{s}$, with remaining parameters
      held at their nominal values.  Here and in
      Figures~\ref{fig:robust-dz} and~\ref{fig:robust-ar} the converged
      iterate (error below $10^{-13}$) is not plotted, so the curve of a
      run with $n$ iterations ends at index $n-1$.}
  ```
- **A1 (only with the short-label Fig. 6), tex 3189–3190:**

  OLD:
  ```
  for seven aspect ratios at fixed
      $\Delta t = 50\,\mathrm{s}$.
  ```
  NEW:
  ```
  for the seven depths
      $L_z \in \{12, 8, 4, 2, 1, 0.5, 0.25\}\,\mathrm{km}$ of
      Table~\ref{tab:mode-damping} at fixed $\Delta t = 50\,\mathrm{s}$.
  ```

## 5. Outside the two sections, noticed while checking

- **tex 2985–2986:** "through Irksome.\nfully discrete system" is missing "The".
- **tex 2986–2994 and the "FGMRES" in captions at 3038, 3048, 3065, 3189, 3291** describe the linear sweeps incorrectly.
  - **Solver:** they ran GMRES with left preconditioning (`LB_time_slice.py:330` with `--richardson`, as launched by `Bash/bash.py:85`).
  - **Stopping test:** residual rtol 1e-8, not a ten-order error reduction.
  - **Error reference:** the converged iterate of an identical warm-up solve (`LB_time_slice.py:463-473`, monitor `:431-438`), not an independent direct solve.
  - Please confirm with the Table 5 agents before editing these.
- **Nonlinear preconditioner (tex 3330):** "preconditioned by $A_\delta^{-1}$" is not exactly what the nonlinear script applies. The missing pressure back-substitution found on 2026-10-16 is still an open decision for Molin.
- **Cosmetic:** panels in the 8/9 and 11/12 pairs sit at different heights (visible on p. 35). Using `\begin{minipage}[t]{0.48\linewidth}` aligns them.

Compute was a few seconds of Python and a PDF page render; no solver runs, no mpiexec. The only file I created outside my folder was the Python cache file `claude_small_patch/__pycache__/restyle_paper_figures.cpython-312.pyc`, which I deleted.

Files are in `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/figure_referee/`:
- `edits.json`
- `edits.diff`
- `ex_article.figure_referee.tex`
- `tex_edits.py`
- `cand/`
- `dump_figure_data.out`
- `replica_check2.out`
- `marker_check.out`
- `legend_overlap.out`
- `legend_overlap_corrected.out`
- `spread_check.out`
