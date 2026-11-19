# Timing set R (Table 2) re-run with the two ranks pinned (mpiexec --bind-to core --map-by socket), 2026-11-18: summary

Runs folder: `/home/mh2520/molin/phd/yr2/data/small_patch/table5_runs_gmres6_R_bound/runs` (24 runs). PCApply = time of the preconditioner applications of the timed solve (stage "Official-Run"; "Second_Run" for set C), maximum over ranks, including the numerical refactorisation. Medians over the repetitions; spread = (max - min)/median.

## Completeness

24 run folders; 24 of the 64 runs of the full script are present. Missing (40): `B_fs_inner1e-10_np2_nx20_nz20`, `B_fs_inner1e-10_np2_nx40_nz40`, `B_mono_inner1e-10_np2_nx20_nz20`, `B_mono_inner1e-10_np2_nx40_nz40`, `C_star_rep1_np2_nx20_nz20`, `C_star_rep2_np2_nx20_nz20`, `C_star_rep3_np2_nx20_nz20`, `C_vanka_rep1_np2_nx20_nz20`, `C_vanka_rep2_np2_nx20_nz20`, `C_vanka_rep3_np2_nx20_nz20`, `T_fs_rep1_np16_nx160_nz160`, `T_fs_rep1_np16_nx20_nz20`, `T_fs_rep1_np16_nx40_nz40`, `T_fs_rep1_np16_nx80_nz160`, `T_fs_rep1_np16_nx80_nz80`, `T_fs_rep2_np16_nx160_nz160`, `T_fs_rep2_np16_nx20_nz20`, `T_fs_rep2_np16_nx40_nz40`, `T_fs_rep2_np16_nx80_nz160`, `T_fs_rep2_np16_nx80_nz80`, `T_fs_rep3_np16_nx160_nz160`, `T_fs_rep3_np16_nx20_nz20`, `T_fs_rep3_np16_nx40_nz40`, `T_fs_rep3_np16_nx80_nz160`, `T_fs_rep3_np16_nx80_nz80`, `T_mono_rep1_np16_nx160_nz160`, `T_mono_rep1_np16_nx20_nz20`, `T_mono_rep1_np16_nx40_nz40`, `T_mono_rep1_np16_nx80_nz160`, `T_mono_rep1_np16_nx80_nz80`, `T_mono_rep2_np16_nx160_nz160`, `T_mono_rep2_np16_nx20_nz20`, `T_mono_rep2_np16_nx40_nz40`, `T_mono_rep2_np16_nx80_nz160`, `T_mono_rep2_np16_nx80_nz80`, `T_mono_rep3_np16_nx160_nz160`, `T_mono_rep3_np16_nx20_nz20`, `T_mono_rep3_np16_nx40_nz40`, `T_mono_rep3_np16_nx80_nz160`, `T_mono_rep3_np16_nx80_nz80`

## Problems

- none

## Set R: tab:rcm (fieldsplit solver, 2 ranks)

| (N_x, N_z) | RCM: its | PCApply n | PCApply s median [min, max] | /N_z | no RCM: its | PCApply s median [min, max] | /N_z |
|---|---|---|---|---|---|---|---|
| (20, 10) | [11] | [11] | 1.56 [1.552, 1.613] | 0.156 | [11] | 2.749 [2.744, 3.043] | 0.275 |
| (20, 20) | [11] | [11] | 2.87 [2.849, 2.885] | 0.143 | [11] | 9.46 [9.427, 9.547] | 0.473 |
| (20, 40) | [10] | [10] | 4.874 [4.854, 4.91] | 0.122 | [10] | 32.08 [31.95, 32.52] | 0.802 |
| (20, 80) | [10] | [10] | 9.273 [9.236, 9.38] | 0.116 | [10] | 127.3 [127.2, 127.7] | 1.59 |

Affine fit of the RCM PCApply against N_z: 0.559 s + 0.1089 s x N_z, R^2 = 0.9992.

Draft rows for tab:rcm (bandwidths unchanged: they depend only on the patch matrices):

```latex
        $(20, 10)$ & 48 & 1.56s & 0.156s & 531 & 2.749s & 0.275s \\
        $(20, 20)$ & 48 & 2.87s & 0.143s & 1051 & 9.46s & 0.473s \\
        $(20, 40)$ & 48 & 4.874s & 0.122s & 2091 & 32.08s & 0.802s \\
        $(20, 80)$ & 48 & 9.273s & 0.116s & 4171 & 127.3s & 1.59s \\
```

## All runs

| tag | rc | wall s | ranks (log) | its | reason | PCApply n | PCApply s | solve s |
|---|---|---|---|---|---|---|---|---|
| R_norcm_rep1_np2_nx20_nz10 | 0 | 19 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 3.043 | 3.059 |
| R_norcm_rep1_np2_nx20_nz20 | 0 | 40 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 9.547 | 9.577 |
| R_norcm_rep1_np2_nx20_nz40 | 0 | 116 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 31.95 | 32.00 |
| R_norcm_rep1_np2_nx20_nz80 | 0 | 423 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 127.3 | 127.5 |
| R_norcm_rep2_np2_nx20_nz10 | 0 | 17 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.749 | 2.765 |
| R_norcm_rep2_np2_nx20_nz20 | 0 | 40 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 9.46 | 9.489 |
| R_norcm_rep2_np2_nx20_nz40 | 0 | 121 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 32.52 | 32.58 |
| R_norcm_rep2_np2_nx20_nz80 | 0 | 424 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 127.7 | 127.8 |
| R_norcm_rep3_np2_nx20_nz10 | 0 | 18 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.744 | 2.76 |
| R_norcm_rep3_np2_nx20_nz20 | 0 | 40 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 9.427 | 9.456 |
| R_norcm_rep3_np2_nx20_nz40 | 0 | 115 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 32.08 | 32.13 |
| R_norcm_rep3_np2_nx20_nz80 | 0 | 422 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 127.2 | 127.3 |
| R_rcm_rep1_np2_nx20_nz10 | 0 | 14 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.613 | 1.627 |
| R_rcm_rep1_np2_nx20_nz20 | 0 | 19 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.885 | 2.914 |
| R_rcm_rep1_np2_nx20_nz40 | 0 | 27 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 4.91 | 4.956 |
| R_rcm_rep1_np2_nx20_nz80 | 0 | 45 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 9.236 | 9.35 |
| R_rcm_rep2_np2_nx20_nz10 | 0 | 14 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.552 | 1.566 |
| R_rcm_rep2_np2_nx20_nz20 | 0 | 19 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.849 | 2.881 |
| R_rcm_rep2_np2_nx20_nz40 | 0 | 27 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 4.854 | 4.901 |
| R_rcm_rep2_np2_nx20_nz80 | 0 | 46 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 9.38 | 9.498 |
| R_rcm_rep3_np2_nx20_nz10 | 0 | 13 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 1.56 | 1.573 |
| R_rcm_rep3_np2_nx20_nz20 | 0 | 18 | 2 | 13/11/11 | CONVERGED_RTOL | 11 | 2.87 | 2.901 |
| R_rcm_rep3_np2_nx20_nz40 | 0 | 27 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 4.874 | 4.92 |
| R_rcm_rep3_np2_nx20_nz80 | 0 | 45 | 2 | 12/10/10 | CONVERGED_RTOL | 10 | 9.273 | 9.388 |
