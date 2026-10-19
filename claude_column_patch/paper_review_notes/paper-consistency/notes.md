# Paper-consistency review of Section 5.7 (sec:single-column), ex_article.tex lines 3829-3950

Read-only review, 2026-09-30. Line numbers refer to ex_article.tex; PDF pages 40-42.

## Must-fix
1. Naming clash. The paper calls the vertex-star patch the "column patch"/"column-patch smoother"
   (lines 195, 539, 1123, 1277-1296, 1341, 1369, 2912, 3017, 3055, 3712-3772, 3969, 3994, 4068,
   4120-4163, 4308-4343) and line 1288 even says "a single column patch highlighted" for the
   two-column star patch. The new subsection uses bare "column patch" for the ONE-column patch at
   lines 3854, 3906, 3938 and "column"/"column/star" in the table header (3921, 3923).
   Fix: write "single-column patch" everywhere in 5.7; header "single column"; change line 1288 to
   "on one star patch highlighted".
2. Smoother KSP. Line 3890 says "six GMRES iterations per level"; line 3841 says "everything else as in
   Appendix B"; Appendix B lines 4358-4360 lists ksp_type:richardson with damping theta, Section 3.4
   lines 1112-1117 and 1185 describe damped Richardson with theta = 0.2. Evidence: committed script
   (claude_column_patch/code_verification/committed_HEAD_Nonlinear_Boussinesq_Irk_SC.py lines 205-208)
   has 'ksp_type': 'gmres', 'ksp_max_it': 6; ksp_view_star_nominal_2ranks.txt lines 165-168 and 225-228
   show "type: gmres, maximum iterations=6" on both levels; LB_time_slice.py --richardson selects
   chebyshev(6); no active ksp_richardson_scale 0.2 in any current script. Correct Appendix B / 3.4, or
   say in 5.7 that the smoother KSP differs.

## Should-fix
3. Line 3897-3901: "the excess behaves as the constant predicts ... grows with N_z (panel (b))". The
   constant is N_z-independent (rho_1 has no N_z; 04_spectra nz sweep, degree 2: 3L lambda_min 0.0660,
   0.0657, 0.0656, 0.0656, 0.0656) and the inner solve alone gives 16, 17, 17, 17. The study lists the
   N_z growth (20, 26, 27, 32, 39 at dt = 50) as "not explained". Reword.
4. Stopping criterion absent (3888-3893). New runs: FGMRES rtol 1e-6, atol 1e-9, snes atol 1e-9, no
   Eisenstat-Walker (study 4.1; test scripts). Appendix A lines 4292-4296 says Section 5.4 uses EW;
   Section 5.2 uses a solution-error criterion. State the criterion; "repeats the sweeps of Section 5.4"
   is also a different shift parametrisation (K delta_thr vs C_1/dt).
5. Line 3941-3945 "Elsewhere it is slower": rows L_z = 12 km (0.93) and N_x = 80 (0.89) at dt = 50 s are
   not slower (rho_1 = 1.1 and 2.0). Say "no faster, and up to four times slower".
6. Line 3942 "rho_1 >~ 10, i.e. dt >~ 5000 s": rho_1 = 10 at dt = 1850 s on the nominal mesh (3.0 at
   1000 s, 72 at 5000 s). Say "dt >~ 2e3 s (5000 s is the first tested value)".
7. Delta x / N_x: paper defines Delta x = L_x/N_x with N_x base-mesh intervals (524-527); rho_1 uses the
   finest of three levels (L_x/(4 N_x)); figure axis says "N_x (coarse mesh)". Say so once at 3866.
8. Line 3865 "star patch stays near 1.3": star two-level 1/lambda_min is 0.8-1.27 (06_k2_constants
   PART 2; 04_spectra 2L 0.788-1.00). Say "stays below 1.3" and name it the two-level constant.
9. Figure caption 3881-3883: 9 of the 10 omitted runs hit the cap; (dt = 1e4 s, K = 1000) is Newton
   divergence. Say "runs that did not converge (the same ten for both patches)".
10. Reproducibility: single-column counts at dt <= 100 s vary +-20% between identical runs (study 4.5;
    N_z = 80 timing repeats 46, 35, 38; the table prints 38). Add a clause.
11. Section 3.4 line 1136 "the choice of a vertex-star patch is essential" and 1145: add "for k = 1;
    for k >= 2 a single-column patch also satisfies the kernel condition with the coarse level, at a
    larger constant (Section 5.7)".
12. Cross-references: roadmap 2900-2901 (add Section 5.7); Discussion 3965-3980 (2 vs up to 20 inner
    FGMRES iterations is the quantitative evidence for "sufficiently accurate approximation");
    Conclusions 4063 "Two implementation findings" (optional third).
13. Table caption 3915 "nominal ... mesh": the paper never names (N_x, N_z) = (20, 20) nominal. State it.
14. delta_thr (3877-3878, 3891-3893): eq (4.63) is two asymptotic inequalities; the code uses
    delta_thr = dt/(2 C_M^2) exactly (claude_paper_shift.py). Define delta_thr once as that number.

## Optional
- K (3877 ff.) collides with the kernel K of eq:kernel-def (1646-1654, 2399).
- "non-zeros" (3844) vs "nonzero" elsewhere (7x); "vertical faces" (3839) vs "vertical facets" (1314).
- "2.5 times fewer" -> give 21 060 against 53 468 (timing logs, rows=364/608 factored matrices).
- "the rest ... is unchanged" (3939-3940): 0.16 vs 0.14 s per iteration -> "nearly unchanged".
- Caption (d) "runs of (a), (b) and the full mesh sweeps" -> "the 170 runs at K = 1 of the aspect-ratio
  and the two mesh sweeps". Underfull vbox on pp. 41-42 from the [t] floats: cosmetic.
- Placement: keep before Discussion (Discussion should synthesise all experiments) and cite it there.

## Verified correct
Labels/bib keys; 18N_z+4 / 30N_z+8, 364/608; LU nnz ratio 2.54; kernel counts (k-1)(kN_z-1),
(2k-1)(kN_z-1), 0 for k=1; 81=81, 9<17; C(rho) 12.28 -> 2.00 (script re-run); 170 runs, star 3-40,
column 10-50; medians 5.20 and 1.11; 32 vs 3; dx 23,24,19,16,15; dz 20..39; K>=100 ratios 0.94-1.06;
inner solve 2 vs <= 20; all ten table rows; 30%/iteration, 2.2x MatSolve, break-even 1.4, gains 18-27%,
factor 4.12; divergence matrix geometry-independent (4e-26, 3e-16); patch composition.
