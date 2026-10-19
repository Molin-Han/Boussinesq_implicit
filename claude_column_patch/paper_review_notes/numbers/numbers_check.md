# Numbers check of Section 5.7 (sec:single-column), ex_article.tex lines 3829-3950

Checked 2026-09-30 against the raw data (recomputed with own scripts from the CSVs; the
two-level constant script re-run in the Firedrake venv). "OK" = matches.

## Patch sizes, kernel dimensions, constants
| statement | value in data | verdict |
|---|---|---|
| 18Nz+4 vs 30Nz+8; 364 vs 608 at Nz=20 | ksp_view_{vanka,star}: rows=364 / 608; patch_dofs report: formula nz(4k^2+k)+2k, nz(8k^2-k)+4k, 58/98 at nz=3 | OK |
| LU factor 2.5 times fewer non-zeros | 53468 / 21060 = 2.54 | OK |
| (k-1)(kNz-1) vs (2k-1)(kNz-1), 0 for k=1 | 01_kernel_dimension_counts: k=2 nz=3: 5/15, nz=5: 9/27, nz=4: 7/21; k=1: 0 / nz-1 | OK |
| 81 = 81 (k=2), 9 < 17 (k=1), 8x3 mesh | 01 and 02 (Firedrake prolongation): dim N_H + column kernels = 81 / 9, dim N_h = 81 / 17 | OK |
| constant 2.0 (rho->inf) to 12.3 (rho->0) | 06_k2_constants: 2.0000, 12.2788; my re-run of my_checks/check_two_level_constant.py: identical | OK |
| star "stays near 1.3" | 1/lambda_min (two-level, k=2): 04_spectra 0.788-1.00 -> 1.00-1.27; 06 PART 2: 0.79-1.24 -> 0.80-1.26 | loose: "below 1.3" |
| rho_1 = m Lz^2/(pi^2 dx^2) | m = 1+dt^2N^2/4 (line 2697); vertical wavenumber pi l/Lz (eq:mode-eigenvalues); dx = finest-level width Lx/(4Nx) in paper_figure.py | OK, but dx level unstated |
| delta_thr = eq (4.63) | claude_paper_shift.shift_threshold = dt/(2C_M^2); dt=50: 1.0966e-8 = 2pi^2 dt/Lx^2; dt=1e5: 1.2337e-7 = 2pi^2/(dt N^2 Lz^2); commands.txt agree | OK |

## Sweeps (summary_iterations.csv, 8 folders; 50/50/70/50 rows)
| statement | recomputed | verdict |
|---|---|---|
| converges wherever star does | non-converged: dt sweep only, 10 runs each, identical set (9 linear cap 100, 1 Newton cap 10 at dt=1e4, K=1000) | OK |
| K=1 totals 10-50 vs 3-40 | 170 K=1 runs (dx,dz,ar): column 10-50, star 3-40; dt sweep K=1: 13-31 / 4-29 | OK |
| 32 vs 3 at Lz=250 m (dt=50) | ar: 3 / 32 | OK (ratio 10.67 = max; absolute excess 29 < 37 of dz Nz=60 dt=100: 13/50) |
| shrinks with Nx, grows with Nz (dt=50) | dx column 23,24,19,16,15 (star 13,13,11,11,11); dz column 20,26,27,32,39 (star 11,13,11,11,11) | OK |
| vanishes for K >= 100 | K=100: ratio 0.94-1.06 (8 runs); K=1000: 1.00 (2); K=10: 1.00-1.67 | OK |
| 170 runs | 50+50+70 = 170 (dt sweep adds 10 more K=1 runs, not in panel d) | OK |
| median 5.2 below 1e-2, 1.1 above 1e2 | bins: 5.20 (13 runs), 3.04 (12), 2.36 (49), 1.52 (27), 1.16 (19), 1.11 (50) | OK |
| MG alone: star 2, column up to 20, bounded | results_sweep_*.csv random rhs: star 2 (1 at K=1000); column 20,17,8,5,3 (K); 17,17,16,9,6 (dt); 11,17,19,19 (Nx); 16,17,17,17 (Nz); 16,17,17,18 (Lz) | OK |

## Table 6 (timing_all_runs.csv, min of 3 repeats of stage Second_Run)
nominal 13/24 4.0543/5.0337 1.2416; dt1000 20/29 5.8222/5.9588 1.0235; dt5000 29/31 8.1154/6.2745 0.773;
dt1e5 18/20 5.3011/4.3527 0.821; thick12000 13/17 3.9792/3.7169 0.934; thin250 3/27 1.3525/5.5745 4.12;
thin250dt1e5 40/44 11.229/8.9397 0.796; nx80 11/14 12.93/11.455 0.886; nz80 11/38 13.178/26.342 1.999;
K100 62/62 16.14/11.806 0.7315.  All ten rows OK. Column counts of the repeats: nz80 46,35,38; thin250 27,30,27.

## Timing text
| statement | recomputed | verdict |
|---|---|---|
| about 30 % per iteration | time/iteration ratio column/star 0.67-0.74 (nominal 0.210/0.312) | OK |
| patch solves 2.2 times cheaper | MatSolve/iteration 0.125 vs 0.056 s (2.24); 2.21-2.27 across configs (MatSolve includes the coarse MUMPS solves, same for both) | OK |
| rest unchanged | (PCApply-MatSolve)/iteration 0.160 vs 0.136 s nominal; 2-15 % lower with column | OK (approx.) |
| within a factor of about 1.4 | 1/0.70 = 1.43; dt1000: 29/20 = 1.45 -> 1.02 | OK |
| rho_1 >~ 10 i.e. dt >~ 5000 s | rho_1(nominal): dt 50: 0.12, 1000: 3.0, 5000: 72, 1e5: 28820; rho_1 = 10 at dt = 1850 s; Lz=250 m dt=1e5: 113 | OK |
| gains 18 to 27 % | 0.82, 0.80, 0.77, 0.73 | OK |
| "Elsewhere it is slower" | Nx=80 (rho_1 = 2.0): 0.89; Lz=12 km (rho_1 = 1.1): 0.93 -> faster, within repeat scatter (up to 27 %) | NOT OK |
| factor of four | 4.12 | OK |

## Figure (Image_patch/column_patch_sweeps.png) vs CSVs
Panel (a): Lz=12 km dt=1: 8/11; Lz=250 m dt=50: 3/32, dt=1e5: 40/44 -- as plotted. (b) values above -- as plotted.
(c) K=100 dt=1000: 125/132; K=1000 dt=10: 51/51; K=0.1 dt=500: 11/55 -- as plotted; omitted runs = the 10 non-converged.
(d) 170 points, max 10.67 -- as plotted. PDF pages 40-42 (Fig. 14, Table 6) match the tex.

## Other
- "six GMRES iterations per level": script and ksp_view (mg_levels_1_: type gmres, maximum iterations=6) -- OK,
  but Appendix B lists "ksp_type": "richardson" and "six Richardson sweeps ... damping factor theta" (pre-existing).
- Section 5 sweeps use delta = C_1/dt, C_1 = 1e-3, dt <= 1e4 s (lines 2993, 3031); the new sweeps use delta = K delta_thr, dt <= 1e5 s.
- sigma (vertical wavenumber here) clashes with sigma_jl (eigenvalues of S) of eq:mode-eigenvalues.
- Not independently recomputed: the star lambda_min tables (04_spectra, 06 PART 2, agent output); timing scatter cause.
