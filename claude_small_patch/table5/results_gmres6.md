# Timing runs with the GMRES(6) level smoother: summary

Runs folder: `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/runs_gmres6` (64 runs). PCApply = time of the preconditioner applications of the timed solve (stage "Official-Run"; "Second_Run" for set C), maximum over ranks, including the numerical refactorisation. Medians over the repetitions; spread = (max - min)/median.

## Completeness

64 run folders; 64 of the 64 runs of the full script are present.

## Problems

- none

## Set R: tab:rcm (fieldsplit solver, 2 ranks)

| (N_x, N_z) | RCM: its | PCApply n | PCApply s median [min, max] | /N_z | no RCM: its | PCApply s median [min, max] | /N_z |
|---|---|---|---|---|---|---|---|
| (20, 10) | [11] | [11] | 1.85 [1.601, 1.919] | 0.185 | [11] | 3.307 [2.779, 3.354] | 0.331 |
| (20, 20) | [11] | [11] | 3.454 [2.875, 3.478] | 0.173 | [11] | 9.483 [9.316, 11.98] | 0.474 |
| (20, 40) | [10] | [10] | 4.902 [4.895, 6.114] | 0.123 | [10] | 42.26 [31.89, 42.66] | 1.06 |
| (20, 80) | [10] | [10] | 9.167 [9.152, 11.78] | 0.115 | [10] | 128.3 [124.7, 172.7] | 1.6 |

Affine fit of the RCM PCApply against N_z: 1.054 s + 0.1010 s x N_z, R^2 = 0.9923.

Draft rows for tab:rcm (bandwidths unchanged: they depend only on the patch matrices):

```latex
        $(20, 10)$ & 48 & 1.85s & 0.185s & 531 & 3.307s & 0.331s \\
        $(20, 20)$ & 48 & 3.454s & 0.173s & 1051 & 9.483s & 0.474s \\
        $(20, 40)$ & 48 & 4.902s & 0.123s & 2091 & 42.26s & 1.06s \\
        $(20, 80)$ & 48 & 9.167s & 0.115s & 4171 & 128.3s & 1.6s \\
```

## Set T: tab:monolithic (16 ranks, FGMRES right, GMRES(6) + star patch + RCM)

| (N_x, N_z) | unknowns | fs its | fs PCApply s median [min, max] (spread) | fs s/unknown | mono its | mono PCApply s median [min, max] (spread) | mono/fs | mono higher by | fs cheaper by |
|---|---|---|---|---|---|---|---|---|---|
| (20, 20) | 32320 | [11] | 1.825 [1.793, 1.836] (2.3%) | 5.65e-05 | [11] | 2.203 [2.083, 2.285] (9.1%) | 1.207 | 20.7% | 17.2% |
| (40, 40) | 128640 | [10] | 5.843 [5.647, 6.117] (8.0%) | 4.54e-05 | [10] | 6.702 [6.645, 6.983] (5.0%) | 1.147 | 14.7% | 12.8% |
| (80, 80) | 513280 | [10] | 23.02 [22.98, 23.14] (0.7%) | 4.48e-05 | [10] | 25.28 [24.98, 25.36] (1.5%) | 1.098 | 9.8% | 9.0% |
| (80, 160) | 1025280 | [9] | 43.01 [42.92, 43.79] (2.0%) | 4.19e-05 | [9] | 45.76 [45.53, 46.4] (1.9%) | 1.064 | 6.4% | 6.0% |
| (160, 160) | 2050560 | [9] | 90.27 [89.63, 92.15] (2.8%) | 4.40e-05 | [9] | 102.0 [101.2, 106.5] (5.2%) | 1.130 | 13.0% | 11.5% |

"mono higher by" (mono/fs - 1) is the comparison of the Cost subsection ("11--17% higher"); "fs cheaper by" (1 - fs/mono) is the one of the Conclusions ("about 10--15% cheaper"). s/unknown divides by the exact number of unknowns (the October text divided by the rounded 3.2e4).

Draft rows for tab:monolithic (unknowns and patch DOFs from the formulas above; check "1.03M"/"2.05M"):

```latex
        $(20, 20)$ & 32k & 608 & 1.82s & 11 & 768 & 2.20s & 11 \\
        $(40, 40)$ & 129k & 1208 & 5.84s & 10 & 1528 & 6.70s & 10 \\
        $(80, 80)$ & 513k & 2408 & 23.02s & 10 & 3048 & 25.28s & 10 \\
        $(80, 160)$ & 1.03M & 4808 & 43.01s & 9 & 6088 & 45.76s & 9 \\
        $(160, 160)$ & 2.05M & 4808 & 90.27s & 9 & 6088 & 102.01s & 9 \\
```

## Set B: the cycle iterated to 1e-10 (2 ranks)

| solver | (N_x, N_z) | outer its (step 1 / warm-up / official) |
|---|---|---|
| fs | (20, 20) | 13/11/11 |
| fs | (40, 40) | 12/10/10 |
| mono | (20, 20) | 13/11/11 |
| mono | (40, 40) | 12/10/10 |

Compare the official counts with the one-cycle counts of set T at (20,20) and (40,40).

## Set C: cost of one outer iteration, star patch against single column (nonlinear solver, 2 ranks)

| patch | runs | applications (= FGMRES its of the timed step) | PCApply s median | s per application median [min, max] |
|---|---|---|---|---|
| star | 3 | [8] | 2.583 | 0.3228 [0.3189, 0.3229] |
| vanka | 3 | [10] | 2.66 | 0.266 [0.2274, 0.2753] |

One application with the single column costs 0.824 times one with the star patch, i.e. it is 18% cheaper (paper and supplement say 'about 30%'). The single column is then faster only where the two counts are within a factor of 1.21.

## All runs

| tag | rc | wall s | ranks (log) | its | reason | PCApply n | PCApply s | solve s |
|---|---|---|---|---|---|---|---|---|
| B_fs_inner1e-10_np2_nx20_nz20 | 0 | 19 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.808 | 2.842 |
| B_fs_inner1e-10_np2_nx40_nz40 | 0 | 56 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 12.19 | 12.32 |
| B_mono_inner1e-10_np2_nx20_nz20 | 0 | 22 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.95 | 3.986 |
| B_mono_inner1e-10_np2_nx40_nz40 | 0 | 60 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 13.93 | 14.06 |
| C_star_rep1_np2_nx20_nz20 | 0 | 161 | 2 | 8 | CONVERGED_ATOL | 8 | 2.551 | 2.828 |
| C_star_rep2_np2_nx20_nz20 | 0 | 20 | 2 | 8 | CONVERGED_ATOL | 8 | 2.583 | 2.867 |
| C_star_rep3_np2_nx20_nz20 | 0 | 21 | 2 | 8 | CONVERGED_ATOL | 8 | 2.583 | 2.865 |
| C_vanka_rep1_np2_nx20_nz20 | 0 | 21 | 2 | 10 | CONVERGED_ATOL | 10 | 2.753 | 3.078 |
| C_vanka_rep2_np2_nx20_nz20 | 0 | 19 | 2 | 10 | CONVERGED_ATOL | 10 | 2.274 | 2.571 |
| C_vanka_rep3_np2_nx20_nz20 | 0 | 21 | 2 | 10 | CONVERGED_ATOL | 10 | 2.66 | 2.969 |
| R_norcm_rep1_np2_nx20_nz10 | 0 | 19 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.779 | 2.794 |
| R_norcm_rep1_np2_nx20_nz20 | 0 | 40 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 9.316 | 9.341 |
| R_norcm_rep1_np2_nx20_nz40 | 0 | 151 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 42.26 | 42.33 |
| R_norcm_rep1_np2_nx20_nz80 | 0 | 427 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 128.3 | 128.5 |
| R_norcm_rep2_np2_nx20_nz10 | 0 | 20 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.307 | 3.325 |
| R_norcm_rep2_np2_nx20_nz20 | 0 | 41 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 9.483 | 9.512 |
| R_norcm_rep2_np2_nx20_nz40 | 0 | 150 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 42.66 | 42.74 |
| R_norcm_rep2_np2_nx20_nz80 | 0 | 573 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 172.7 | 172.8 |
| R_norcm_rep3_np2_nx20_nz10 | 0 | 22 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.354 | 3.371 |
| R_norcm_rep3_np2_nx20_nz20 | 0 | 50 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 11.98 | 12.02 |
| R_norcm_rep3_np2_nx20_nz40 | 0 | 116 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 31.89 | 31.95 |
| R_norcm_rep3_np2_nx20_nz80 | 0 | 417 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 124.7 | 124.8 |
| R_rcm_rep1_np2_nx20_nz10 | 0 | 118 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.919 | 1.938 |
| R_rcm_rep1_np2_nx20_nz20 | 0 | 20 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.875 | 2.904 |
| R_rcm_rep1_np2_nx20_nz40 | 0 | 33 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 6.114 | 6.18 |
| R_rcm_rep1_np2_nx20_nz80 | 0 | 46 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 9.152 | 9.256 |
| R_rcm_rep2_np2_nx20_nz10 | 0 | 16 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.601 | 1.616 |
| R_rcm_rep2_np2_nx20_nz20 | 0 | 22 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.454 | 3.489 |
| R_rcm_rep2_np2_nx20_nz40 | 0 | 29 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 4.902 | 4.96 |
| R_rcm_rep2_np2_nx20_nz80 | 0 | 46 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 9.167 | 9.259 |
| R_rcm_rep3_np2_nx20_nz10 | 0 | 16 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.85 | 1.868 |
| R_rcm_rep3_np2_nx20_nz20 | 0 | 22 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.478 | 3.514 |
| R_rcm_rep3_np2_nx20_nz40 | 0 | 28 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 4.895 | 4.953 |
| R_rcm_rep3_np2_nx20_nz80 | 0 | 55 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 11.78 | 11.91 |
| T_fs_rep1_np16_nx160_nz160 | 0 | 338 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 92.15 | 93.02 |
| T_fs_rep1_np16_nx20_nz20 | 0 | 20 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 1.825 | 1.845 |
| T_fs_rep1_np16_nx40_nz40 | 0 | 34 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 6.117 | 6.186 |
| T_fs_rep1_np16_nx80_nz160 | 0 | 171 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 43.79 | 44.22 |
| T_fs_rep1_np16_nx80_nz80 | 0 | 94 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 23.14 | 23.37 |
| T_fs_rep2_np16_nx160_nz160 | 0 | 325 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 89.63 | 90.5 |
| T_fs_rep2_np16_nx20_nz20 | 0 | 20 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 1.836 | 1.856 |
| T_fs_rep2_np16_nx40_nz40 | 0 | 32 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 5.647 | 5.707 |
| T_fs_rep2_np16_nx80_nz160 | 0 | 169 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 42.92 | 43.37 |
| T_fs_rep2_np16_nx80_nz80 | 0 | 92 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 22.98 | 23.21 |
| T_fs_rep3_np16_nx160_nz160 | 0 | 328 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 90.27 | 91.15 |
| T_fs_rep3_np16_nx20_nz20 | 0 | 19 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 1.793 | 1.813 |
| T_fs_rep3_np16_nx40_nz40 | 0 | 33 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 5.843 | 5.906 |
| T_fs_rep3_np16_nx80_nz160 | 0 | 167 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 43.01 | 43.44 |
| T_fs_rep3_np16_nx80_nz80 | 0 | 92 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 23.02 | 23.25 |
| T_mono_rep1_np16_nx160_nz160 | 0 | 369 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 106.5 | 107.4 |
| T_mono_rep1_np16_nx20_nz20 | 0 | 19 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 2.285 | 2.308 |
| T_mono_rep1_np16_nx40_nz40 | 0 | 35 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 6.983 | 7.046 |
| T_mono_rep1_np16_nx80_nz160 | 0 | 179 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 46.4 | 46.85 |
| T_mono_rep1_np16_nx80_nz80 | 0 | 98 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 25.28 | 25.52 |
| T_mono_rep2_np16_nx160_nz160 | 0 | 355 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 101.2 | 102.0 |
| T_mono_rep2_np16_nx20_nz20 | 0 | 18 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 2.083 | 2.108 |
| T_mono_rep2_np16_nx40_nz40 | 0 | 34 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 6.702 | 6.767 |
| T_mono_rep2_np16_nx80_nz160 | 0 | 173 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 45.53 | 45.95 |
| T_mono_rep2_np16_nx80_nz80 | 0 | 98 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 25.36 | 25.6 |
| T_mono_rep3_np16_nx160_nz160 | 0 | 360 | 16 | 11/9/9 | CONVERGED_RTOL | 9 | 102.0 | 102.9 |
| T_mono_rep3_np16_nx20_nz20 | 0 | 19 | 16 | 13/11/11 | CONVERGED_RTOL | 11 | 2.203 | 2.225 |
| T_mono_rep3_np16_nx40_nz40 | 0 | 34 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 6.645 | 6.708 |
| T_mono_rep3_np16_nx80_nz160 | 0 | 175 | 16 | 12/9/9 | CONVERGED_RTOL | 9 | 45.76 | 46.19 |
| T_mono_rep3_np16_nx80_nz80 | 0 | 97 | 16 | 12/10/10 | CONVERGED_RTOL | 10 | 24.98 | 25.22 |
