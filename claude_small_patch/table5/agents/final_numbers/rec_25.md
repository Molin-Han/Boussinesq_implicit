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

| (N_x,N_z) | FS left | FS right | MONO (as run) left | MONO (as run) right | MONO (as run), GMRES(6) smoother, FGMRES right |
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

| (N_x,N_z) | MONO RCM + PETSc LU | MONO natural + PETSc LU | MONO MUMPS | FS RCM + PETSc LU | FS MUMPS |
|---|---|---|---|---|---|
| (20,20) | 7/6/6 | 7/6/6 | 7/6/6 | 11/13/13 | 11/13/13 |
| (40,40) | 7/7/7 | 7/7/7 | 7/7/7 | 12/13/13 | 12/13/13 |

(Timings of Set Y are not comparable: it ran next to other jobs. MUMPS patch solves are 4-6 x slower than the banded PETSc LU with RCM.)

#### (e) Larger meshes, 2026-06-01 configuration (Set C, GMRES left)

| (N_x,N_z) | FS 2 ranks | FS 16 ranks | MONO (as run) 2 ranks | MONO (as run) 16 ranks | FS PCApply 16 ranks | MONO PCApply 16 ranks |
|---|---|---|---|---|---|---|
| (80,80) | 12/13/13 | 12/13/13 | 17/100/100 (cap) | 17/100/100 (cap) | 28.86 | 218.76 |
| (80,160) | - | 12/13/13 | - | - | 70.74 | - |

The rest of Set C (FS and MONO at (80,160) on 2 ranks, MONO at (80,160) on 16 ranks, both solvers at (80,80) in serial) was dropped once the cause was known: each MONO run hits the cap and takes 40-100 min; the serial question is answered by (a) and (f).

#### (f) The corrected monolithic solver (Set F; `-firedrake_2_monoshiftpc_pc_use_amat false`)

| (N_x,N_z) | FS left | FS (pmat) left | MONO fixed left | FS right | MONO fixed right | MONO fixed serial | MONO fixed 16 ranks | MONO as run left |
|---|---|---|---|---|---|---|---|---|
| (10,10) | 11/12/12 | 11/12/12 | 10/12/12 | 13/11/11 | 12/11/11 | - | - | 5/5/5 |
| (20,20) | 11/13/13 | 11/13/13 | 11/13/13 | 13/11/11 | 13/11/11 | 11/13/13 | 11/13/13 | 7/6/6 |
| (40,40) | 12/13/13 | 12/13/13 | 12/13/13 | 13/11/11 | 13/11/11 | 12/13/13 | 12/13/13 | 7/7/7 |
| (40,80) | 12/13/13 | 12/13/13 | 12/13/13 | 13/11/11 | 13/10/10 | - | 12/13/13 | 7/7/7 |

| (N_x,N_z) | FS with exact LU on S (serial, GMRES left) | FS inner 1e-10 (FGMRES right) | MONO fixed, GMRES(6) smoother, FGMRES right | MONO as run, same |
|---|---|---|---|---|
| (20,20) | 11/12/12 | 13/11/11 | 13/11/11 | 2/2/2 |
| (40,40) | 11/12/12 | 12/10/10 | 12/10/10 | 2/2/2 |

#### (h) The 2026-06-01 monolithic configuration on the new meshes (Set H, 16 ranks)

| (N_x,N_z) | MONO as run, FGMRES right | MONO fixed, FGMRES right | FS, FGMRES right |
|---|---|---|---|
| (20,20) | 7/6/6 | 13/11/11 | 13/11/11 |
| (40,40) | 8/8/8 | 13/11/11 | 13/11/11 |
| (80,80) | 100/100/100 (cap) | 12/10/10 | 13/11/11 |

MONO as run with MUMPS patch LU at (80,80), GMRES left, 16 ranks: -; with a residual monitor (FGMRES right): 100 (cap).

