# paper_text agent report (workflow wf_68296924-398, 2026-10-19)

Paper-text review of the uncommitted ex_article.tex diff, most severe first. Nothing outside the scratch directory was changed. All 20 replacements below were applied together to a copy and built with `latexmk -norc`: 0 errors, 0 undefined references, overfull boxes go from 13 to 12, and Table 5 stays on p.41 with its discussion. The exact old/new pairs are in `edits.json` and the full diff in `proposed_edits.diff`.

**1. Wrong percentage, Table-5 paragraph, l.3810–3814.** "the cumulative PCApply time of the fieldsplit Schur approach is lower by 11–17%" uses the wrong base.
- Evidence: FS is lower than mono by 14.6, 13.5, 9.6 and 11.3% (so 10–15%). Mono is higher than FS by 17.1, 15.6, 10.7 and 12.7% (so 11–17%).
- Fix: see issue 4. Its replacement sentence reads "its cumulative PCApply time is $11$--$17\%$ higher than that of the fieldsplit Schur approach at every resolution except $(80, 80)$, where the monolithic cycle saves an iteration and the difference is $1\%$."

**2. Wrong fact, Table 5 caption, l.3832–3834.** "...RCM-reordered patch factorisations on all three levels" is wrong. The coarsest level is a direct LU (MUMPS); Chebyshev(6) with patch smoothing runs on levels 1 and 2 only.
- Evidence: T_*/slice3D.txt shows "levels=3", "mg_coarse ... lu ... mumps", and Chebyshev on levels 1 and 2.
- Fix, whole caption tail from "cycle per application in both (" to "medians of three runs.}":
```
    cycle per application in both (three levels: Chebyshev(6) smoothing
    with the column-star patch and RCM-reordered patch factorisations on
    the two finer levels, a direct solve on the coarsest).  FGMRES with
    right preconditioning to a relative residual of $10^{-8}$; aspect
    ratio $\varepsilon = 1.3\times10^{-2}$, $\Delta t = 100\,\mathrm{s}$,
    $\delta = 10^{-4}$ ($C_1 = \delta\,\Delta t = 10^{-2}$), $16$ MPI
    processes.  $N_x$ is the number of columns of the coarsest multigrid
    level (the finest has $4N_x$); the patches have the same size on
    both smoothed levels; times are medians of three runs.}
```
- This also fixes two smaller points. "patch DOFs are those of the finest level" implies coarser levels differ; both smoothed levels have rows=3048 (mono) and 2408 (FS) at (80,80). And δ=1e-4 at Δt=100 s means C₁=1e-2, while Sec. 5.2 says δ=C₁/Δt "throughout" with nominal C₁=1e-3.

**3. Wrong fact in the edited sentence, RCM paragraph l.3737–3739, plus Table 4 caption l.3756.** "the wall-clock cost of a single preconditioner application (PCApply)" is wrong. The Table 4 numbers are the Official-Run PCApply total over 14 applications (13 GMRES-left iterations plus 1), and that total includes the patch refactorisation.
- Evidence: data_firedrake26_reordering, all 8 logs show "with 2 processes" and "PCApply 14". For example 2.6279 s at (20,20) with RCM.
- Fix, l.3737–3739:
```
Table~\ref{tab:rcm} reports the patch matrix bandwidth and the
wall-clock time spent in the preconditioner (PCApply) during one GMRES
solve, over its $14$ applications, with and without RCM reordering,
for a fixed aspect ratio
```
- Fix, Table 4 caption, from the caption start up to "Aspect ratio $\varepsilon = 1.3\times10^{-2}$":
```
    \caption{Effect of RCM reordering on the column-patch LU
    factorisations. PCApply is the wall-clock time of the $14$
    preconditioner applications of one GMRES solve on $2$ MPI
    processes, including the numerical refactorisation of the patch
    matrices. Aspect ratio $\varepsilon = 1.3\times10^{-2}$
```
- Fix, l.3748, stale header name: `$\mathrm{Time(PCApply)}$` → `PCApply`.

**4. Overclaim, l.3807–3808.** "The cost is therefore decided by the patches" is not what the profile shows.
- Evidence: patch-solve time per iteration (medians of 3 runs, Official-Run stage) is only −2 to +10% higher for mono. MatSolve: +2.5, +10.3, +5.8, +2.2, +6.0%. PCApplyOnBlocks: −1.5, +9.0, +7.0, +4.8, +5.4%. Total PCApply per iteration is +11 to +17%.
- Fix, replace l.3802 "therefore gives exactly..." through l.3814 "...difference is $1\%$.":
```
therefore gives exactly the fieldsplit patch problem on $S$; the
coarse-grid corrections agree in the velocity, because the spaces are
nested and the cells affine, so that the rediscretised coarse operators
are the Galerkin ones and the divergence commutes with the
prolongation; and, since PETSc's eigenvalue estimates give both
smoothers the same Chebyshev interval ($[0.2, 2.2]$ to four digits),
one monolithic cycle on $A_\delta$ returns the same velocity as one
fieldsplit cycle.  The two differ only in the pressure, which the
Schur factorisation recovers exactly from the shifted divergence
equation~\eqref{eq:shifted-div} and the monolithic cycle returns with
a small residual.  The difference in cost therefore comes from the
larger systems the monolithic cycle works with: its patches also carry
the pressure unknowns of their two columns, $38N_z + 8$ against
$30N_z + 8$ unknowns, and its smoother, residuals and grid transfers
act on the pressure as well.  An application of the monolithic cycle
is $11$--$17\%$ more expensive, and its cumulative PCApply time is
$11$--$17\%$ higher than that of the fieldsplit Schur approach at every
resolution except $(80, 80)$, where the monolithic cycle saves an
iteration and the difference is $1\%$.
```
- This also makes the coarse-grid argument explicit (Galerkin coarse operators, commuting divergence). It states the Chebyshev-interval condition, which holds: [0.2, 2.2] in every T run, and 2.20004 vs 2.20001 on one level at (160,160). And it uses the paper's own term "shifted divergence equation" (eq:shifted-div) instead of "shifted continuity equation".

**5. Overclaim, conclusion l.4066–4073.** "converges in the same number of iterations as the Schur reduction, as it must" is too strong.
- Evidence: the counts differ at (80,80) (11 vs 10). Equal velocity output does not force equal counts, because the pressure differs; the theory agent's 8×8 run at δ=1e-3 gives 22 vs 24. The patch identity alone is also not the whole reason.
- The three commented lines that follow ("% so that the cumulative ... lower for the monolithic solver ...") contradict the new text and should be deleted.
- Fix, replace l.4066–4073 including the comment lines:
```
monolithic column-patch smoother applied directly to $A_\delta$
returns the same velocity as the Schur reduction in each application,
since eliminating the pressure from its patches gives the fieldsplit
patches, and converges in the same number of iterations to within one;
the fieldsplit patches are about $20\%$ smaller, and the Schur
reduction is the cheaper of the two, by $10$--$15\%$ in preconditioner
time at four of the five resolutions tested.
```

**6. Inconsistency, Appendix B l.4354–4360 (and Sec. 3.4 l.1186).** The appendix shows a Richardson(6) smoother and says this is "the sense in which the comparison of Section 5.6 ... is untuned on both sides". The new Table 5 caption says Chebyshev(6), which is what ran (`--richardson` sets the Chebyshev type; ksp_view confirms).
- Fix, from "with the damping factor $\theta$.  Neither has been" to "...solved directly.":
```
with the damping factor $\theta$.  The timings of
Section~\ref{sec:timing} use \opt{ksp_type:chebyshev} with
\opt{ksp_max_it:6} instead, six Chebyshev steps on the interval
estimated by PETSc, for the fieldsplit and the monolithic smoother
alike.  Neither choice has been tuned.  The coarse grid is solved
directly.
```
- Optional, l.1186–1187: replace "factor of $\theta = 0.2$. " with `factor of $\theta = 0.2$ (the timings of Section~\ref{sec:timing} use a Chebyshev iteration instead; Appendix~\ref{app:petsc-options}). `
- If the Sec. 5.3 sweeps also used `--richardson` (the claude_paper_linear_* commands do), say "Sections 5.3 and 5.6". I could not check the commands behind the original data_linear_*_deg2 runs.

**7. Notation conflict for N_x.** The new caption defines N_x as the column count of the coarsest level. Sec. 2.2 (l.524–526) defines Δx = L_x/N_x, and Sec. 3.4 (l.1095, 1108) sums over N_x base vertices.
- Evidence: in the code `--nx` is the coarse mesh and `--refinement` defaults to 2, and no Section-5 command overrides it.
- Fix: append after l.3024 "...for better ordering of the sparse matrix.":
```
The multigrid hierarchy has three levels: a base mesh of $N_x$
intervals, refined twice in the horizontal direction, with $N_z$ layers
on every level.  The mesh sizes $(N_x, N_z)$ quoted in this section
refer to the coarsest level, so that the finest mesh, on which the
system is solved, has $4N_x \times N_z$ cells and
$\Delta x = L_x/(4N_x)$.
```

**8. PCApply includes the setup, l.3792–3794.** In the T runs, PCSetUp (numerical refactorisation of the patch and coarse LUs) sits inside PCApply and is 10–23% of it.
- Evidence: KSPSolve − PCApply ≤ 1.0 s, while PCSetUp is 20.8 s at FS (160,160). Without the setup, mono is still 11–16% more per application, so the conclusions hold.
- Fix, l.3792–3798 from "degrees of freedom per patch on the finest level," to "This is to be expected.":
```
degrees of freedom per patch, the wall-clock time spent in the
preconditioner (PCApply) during one FGMRES solve, which includes the
numerical refactorisation of the patch and coarse-grid matrices, and
the iteration count to a relative residual of
$10^{-8}$, for both approaches with one multigrid cycle per application
and the same smoother.  The two converge in the same number of
iterations, except at $(80, 80)$, where the monolithic cycle needs one
fewer.  That they agree is to be expected.
```
- "This is to be expected" read as if it referred to the exception; the new wording removes that.

**9. Second paragraph, l.3816–3826.** Four problems:
- It should use the paper's notation A_0.
- "it approximates the inverse of the unshifted operator" is supported only on 16×4 (2.4% from A_0⁻¹r), and the cycle fails at (80,80).
- "coarser meshes" can be confused with multigrid levels.
- "is handed" is informal.
- Fix, whole paragraph:
```
A monolithic cycle in which only the patch factorisations and the
coarser levels use $A_\delta$, while the finest-level smoother and
residual use the unshifted operator $A_0$, is a different
preconditioner.  It is what PETSc's multigrid does by default when it
is given $A_0$ as the operator and $A_\delta$ only as the
preconditioning matrix, as it is inside Firedrake's auxiliary-operator
preconditioner; the option \opt{pc_use_amat:false} makes it use
$A_\delta$ throughout, as in Table~\ref{tab:monolithic}.  That cycle
is a multigrid method for $A_0$ rather than an approximation of
$A_\delta^{-1}$.  It needs fewer iterations on the two smaller meshes
($6$ and $8$ at $(20, 20)$ and $(40, 40)$), fewer even than the exact
inverse of $A_\delta$ ($11$ and $10$), but it does not converge within
$100$ iterations at $(80, 80)$, so it is not robust.
```

**10. Wrong citation, existing text next to the edits, l.274 and l.3785–3787.** farrellMitchellWechsung2019 is an augmented-Lagrangian block (Schur) preconditioner with star-patch multigrid on the augmented velocity block. That is the analogue of the fieldsplit approach, not of the monolithic one, and it is not discussed in Sec. 3.4. This is from my knowledge of that paper; please check it.
- Fix, l.3785–3787:
```
monolithic construction is a column-patch instance of the monolithic
patch relaxation of the PCPATCH framework~\cite{PCPATCH} of
Section~\ref{sec:inner-solvers}.
```
- Fix, l.274: `~\cite{farrellMitchellWechsung2019}` → `~\cite{PCPATCH}`. FMW2019 then drops out of the bibliography unless it is re-cited as the augmented-Lagrangian analogue.

**11. Introduction overstates the cost, l.275.** "carry a high cost per application because of their size" does not match the measured 11–17%.
- Fix: `Monolithic patches are robust but larger, and therefore more expensive to apply. `

**12. Minor.**
- l.3791: "$2.1\times10^{6}$" does not match the table's 2.05M. Fix: `$3.2\times10^{4}$ to $2.05\times10^{6}$ unknowns`.
- l.1315–1316, Sec. 3.4: "horizontal-velocity degrees of freedom on their vertical facets" should read "on the vertical facet they share". The new elimination argument needs exactly this, and the 30N_z+8 count confirms only the shared facet is in the patch.

**13. Optional, existing text.**
- l.3745–3747: the per-layer cost "is essentially constant" while the table shows it falling from 0.169 to 0.100. Replace from "bandwidth is independent of $N_z$ and the PCApply cost per layer," up to "(a least-squares fit of PCApply against $N_z$" with:
```
bandwidth is independent of $N_z$ and the PCApply cost grows
linearly in $N_z$: the cost per layer, reported in the right-most
column of each block, falls from $0.169\,\mathrm{s}$ to
$0.100\,\mathrm{s}$ as a fixed overhead is amortised
(a least-squares fit of PCApply against $N_z$
```
  The R²=0.998 fit is correct.
- l.4065: "the cost per layer is constant" → "the cost grows only linearly".
- l.3735: "all FGMRES iterations" → "all Krylov iterations" (the Table 4 runs used GMRES).
- l.3021–3022: replace `\texttt{ASMStarPC} (\texttt{pc\_star\_construct\_dim:0})` with `\opt{ASMStarPC} (\opt{pc_star_construct_dim:0})`. This removes the 55.5pt overfull box.

**14. Layout.** If the text grows, as it does with the fixes above, Table 5 floats to p.42 inside Sec. 5.7. Moving the `\begin{table}...\end{table}` block to just before "Table~\ref{tab:monolithic} reports, ..." keeps it on p.41 (tested).

**Checks that passed**
- **Maths.** The velocity equivalence is correct; I re-derived it. In the variables (e_v, η) with η = e_p + δ⁻¹M_p⁻¹Gᵀe_v, the monolithic patch correction is (S_i⁻¹R_iS e_v, R_iᵖη). The smoother is diag(B_S S, 2I); the 2I part needs BASIC ASM, which asm.py sets, but the velocity part holds for any ASM variant. The coarse correction is diag(I−P S_H⁻¹PᵀS, I−Π_H), and this carries through the FMG cycle.
- **Refutation attempts that failed.** The constant-pressure null space attached to the monolithic auxiliary operator does not change the velocity, because a constant pressure shift leaves the reduced residual unchanged. Dirichlet rows are consistent (the audit measured 1e-16). The Chebyshev intervals match. Equal counts to within one follows from a pressure-only perturbation of about 1e-3.
- **Table 5 numbers.** Counts, median PCApply, patch DOFs (38N_z+8 and 30N_z+8), unknowns, per-application 11–17%, 1% at (80,80), "about 20% smaller" (21%), 16 ranks, FGMRES right with rtol 1e-8 on the unpreconditioned norm, FMG with 3 levels, and the second paragraph's 6/8, 11/10 and the 100-iteration cap all match the run logs.
- **Timing introduction and Table 4.** 2 and 16 processes, ε=1.3e-2, L_x=300 km, L_z=4 km, Δt=100 s and δ=1e-4 all match the logs.
- **LaTeX on the current diff.** No errors, no undefined references or citations, no overfull box in the edited region; the 13 overfull boxes elsewhere were there before. Table 4 now fits. Table 5 renders within the text width on p.41.

**Side note for the record (not the paper):** Sections 0 and 2.7 of claude_paper_update_and_table5.md also treat PCApply as pure application time, although it includes the refactorisation (issue 8). Section 2.8 item 6 says the Appendix "still holds"; issue 6 shows it does not.

Files are in /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_paper_text/:
- edits.json
- proposed_edits.diff
- scripts/apply_proposed_v2.py
- build/proposed/ex_article.pdf
- build/page-40.png
- build/page-41.png
- build/proposed_page-40.png
- build/proposed_page-41.png
