# numbers agent report (workflow wf_68296924-398, 2026-10-19)

Verification of Sections 0 and 2 of /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_paper_update_and_table5.md and the Table-5 edits in `git diff HEAD` of ex_article.tex. Every count and time was re-derived from runs/*/run.log, log.txt and slice3D.txt.

**ISSUES (most severe first)**

1. **[HIGH] The caption and record say the patch smoothing is on all three levels; the coarsest level is a direct LU.**
   - Location: paper ex_article.tex:3831-3834 (caption: "Chebyshev(6) smoothing with the column-star patch and RCM-reordered patch factorisations on all three levels"); record line 397, Section 2.7.
   - Evidence:
     - LB_time_slice.py:271-274 and :377-380 set `mg_coarse` to preonly + lu.
     - Every T-run slice3D.txt shows "Coarse grid solver -- level 0 … type: preonly … type: lu … type: mumps". Patch sizes are printed for 2 levels only.
     - The paper's own Appendix (ex_article.tex:4359) says "The coarse grid is solved directly".
   - Fix (caption): "(three-level full multigrid: Chebyshev(6) smoothing with the column-star patch and RCM-reordered patch factorisations on the two finer levels, direct LU on the coarsest)". Make the same change at record line 397.

2. **[HIGH] The cumulative-cost percentage is quoted the wrong way round.**
   - Location: paper ex_article.tex:3811-3812 ("…the cumulative PCApply time of the fieldsplit Schur approach is lower by 11–17%"); record line 45 (Section 0), line 425 (Section 2.7), line 445 (Section 2.8 item 4).
   - Evidence from the raw medians:
     - 1 − fs/mono = 14.3, 13.4, 9.6, 11.3 % at (20,20), (40,40), (80,160), (160,160). So the fieldsplit time is 10–14 % lower.
     - 11–17 % (16.6, 15.5, 10.7, 12.7 %) is mono/fs − 1, i.e. how much more the monolithic run costs.
   - Fix (paper): "…an application of the monolithic cycle is 11–17% more expensive, and so is its cumulative PCApply time at every resolution except (80,80), where … the difference is 1%". Alternatively keep the sentence and write "lower by 10–14%". Apply the same fix in the record.

3. **[HIGH] Section 0 claims a Table-5 commit that does not exist.**
   - Location: record line 47, "Commit COMMIT2 (pushed)".
   - Evidence: "COMMIT2" is an unfilled placeholder. `git status` shows `M ex_article.tex`, and HEAD = origin/main = 0778f01.
   - Fix: commit and push, then insert the real hash; or say the edits are uncommitted.

4. **[MEDIUM] The paper now contradicts itself on the smoother.**
   - Location: Appendix ex_article.tex:4354-4359 says the Section 5.6 comparison uses `ksp_type:richardson`, `ksp_max_it:6` with damping θ, "untuned on both sides". The new Table 5 caption says Chebyshev(6), which is what actually ran (LB_time_slice.py:284-286; ksp_view shows "type: chebyshev").
   - The record's 2.8 item 6 (line 458, "remark still holds") hides this.
   - Fix: describe Chebyshev(6) with PETSc's eigenvalue estimate in the Appendix, or list it as an open inconsistency in 2.8.

5. **[MEDIUM] The claim that the Chebyshev estimate is partition-independent and exactly 2 is false for the old monolithic run.**
   - Location: record lines 365-366 (Section 2.6(a)); lines 198-199 (Section 2.2 item 5, "[0.2, 2.2] on every level"); line 256 (Section 2.4 multigrid row).
   - Evidence:
     - MONO as run, finest level, largest-eigenvalue estimate at (20,20): 2.00776 (1 rank), 2.00632 (2 ranks), 2.00411 (16 ranks) (runs/A_mono_np{1,2,16}_nx20_nz20/slice3D.txt). It is 2.001–2.018 in every run of that configuration.
     - The minimum estimates differ with the partition for every solver.
     - The "noisy" right-hand side is a hash of the PETSc global index (petsc/include/petsc/private/kspimpl.h:260-280), and that index depends on the partition.
     - At (160,160) the finest-level estimate is 2.00003 (fs, fspmat) and 2.00001 (monofix).
   - Fix: "For FS and the corrected MONO the estimate is 2 (to ≤3e-5) on every level and partition. For MONO as run the finest-level estimate depends on the partition (2.0078/2.0063/2.0041 on 1/2/16 ranks), yet the counts are identical."

6. **[MEDIUM] "Does not converge at all at (80,80)" overstates the failure.**
   - Location: record line 35 (Section 0). Also line 515 (Section 4, outside my scope): "3 x 100 iterations with left".
   - Evidence:
     - GMRES left, the 2026-06-01 configuration: 17/100/100 on 2 and 16 ranks (C_mono_*). The first solve converges.
     - FGMRES right: 100/100/100 (H_mono_right).
     - The finished monitor run H_mono_monitor_right_np16_nx80_nz80 reaches relative residual 9.0e-8 at iteration 17 of solve 1, then 2.0e-8 after 100; solves 2 and 3 end at 4.7e-7.
   - Fix: "stalls at (80,80): 17/100/100 (GMRES left), 100/100/100 (FGMRES right; residual reduced only 2e-8 and 4.7e-7 after 100 iterations)".

7. **[MEDIUM] The machine description is wrong and the "one run at a time" claim is not quite true.**
   - Location: record line 260 (Section 2.4 machine row); line 516 (Section 4, "one run at a time").
   - Core count: lscpu reports 2 sockets × 8 cores × 2 threads, i.e. 16 physical cores and 32 hardware threads, not 32 cores.
   - Overlapping jobs:
     - A 2-rank nonlinear job (claude_small_patch/diag_item3/runs/G1_offdiag_amat_dt5000.0, 09:23:40–09:25:21) overlapped the official stage of T_fs_rep3_np16_nx160_nz160 (finished 09:24:34). That run's PCApply of 95.30 s is the maximum of the three; the median is unaffected.
     - C_fs_np16_nx80_nz160 (07:23:01–07:27:05; table (e) value 70.74 s) overlapped Y_mono_rcm_mumps_np2_nx80_nz80 and F_monofix_np2_nx10_nz10. Its cost is 5.05 s per application, against 4.46 s for T_fs on the same mesh.
   - Fix: say "16-core (32-thread)", and note both overlaps.

8. **[MEDIUM-LOW] The (80,80) difference is resolved, contrary to the record.**
   - Location: record line 431 ("the 1 % at (80,80) is not [larger than the scatter]") and line 399 ("fspmat … timing within the run-to-run scatter").
   - Evidence at (80,80):
     - The scatter is 0.2–0.7 %, and the ranges do not overlap: fs [23.32, 23.38], monofix [23.52, 23.69]. So the 1 % difference is real.
     - fspmat [23.82, 23.92] is 2.1 % above fs, outside that mesh's scatter.
   - Fix: reword both sentences accordingly.

9. **[LOW-MEDIUM] The experiment register (Section 2.4) is incomplete or imprecise.**
   - "MONO right" in Set A′ and "MONO fixed right" in Set F ran FGMRES (the monolithic default). "FS right" (Ap_fs_right) ran GMRES with right preconditioning. ksp_view confirms both, but tables (b) and (f) just say "right".
   - Runs missing from the register:
     - Y_mono_natural_right_np2_nx{20,40}: 7/6/6 and 8/8/8.
     - Y_mono_rcm_mumps_np2_nx80_nz80: stopped by hand at 07:31.
     - H_mono_monitor_right_np16_nx80_nz80.
   - The "agents" row omits several theory_blind runs (theory_blind/runs/summary_iterations.csv):
     - shifts δ = 1e-3, 1e-5, 1e-6 at (8,8) (runs E6);
     - meshes (4,8) and (8,4) (runs E7);
     - runs with the Chebyshev interval fixed at [0.2, 2.2].
   - The E6 runs show FS against the corrected MONO at 17/22/22 vs 15/24/24 (δ = 1e-3) and 5/5/5 vs 4/4/4 (δ = 1e-6), both GMRES left. This bears on Section 4's "(to within one) the same counts".

10. **[LOW] Wrong cross-references.**
    - Line 230: "(Section 2.5 (f))" should be (e).
    - Line 202: "Section 2.7" contains no agent checks; they are in the 2.4 "agents" row and Section 3.
    - Section 2.5 jumps from (f) to (h).

11. **[LOW] The 1.7e-23 figure is misattributed.**
    - Location: line 220, "equal to the auxiliary operator to 1.7e-23".
    - Evidence: numerics agent's 1.68e-23 is `fine_matschur_vs_dense_schur_relF`, i.e. MatSchurComplement against the dense Schur complement. The auxiliary operator against the dense Schur complement is 1.03e-16.
    - Fix: say what was compared, and quote 1.0e-16 for the auxiliary operator.

12. **[LOW] Wrong patch-LU accuracy figure.**
    - Location: line 379, "audit measured 1.3e-15".
    - Evidence: the monolithic patch LU is 1.1e-15 on both levels (code_audit/audit/analyse_nx4.txt). 1.2e-15 and 1.3e-15 are the fieldsplit patches.

13. **[LOW] The MUMPS slowdown range is too narrow.**
    - Location: line 326, "MUMPS … 4-6 x slower".
    - Evidence: PCApply ratios are MONO 3.9× / 6.8× and FS 2.7× / 6.1× at (20,20) / (40,40). Use "3–7×".

14. **[LOW] Run times for the dropped Set C runs are overstated.**
    - Location: line 335, "each MONO run hits the cap and takes 40-100 min".
    - Evidence: the cap runs that executed took 497 s (16 ranks) and 1192 s (2 ranks). The 100 min figure is the MUMPS 2-rank estimate in run_table5_pivot.log.

15. **[LOW] Section 0's list of old-table errors is incomplete.**
    - It omits that the fieldsplit counts at (40,40) and (40,80) were 13, printed as 12. Section 2.1 has this.

16. **[LOW] "Serial … for every solver and mesh" is too broad.**
    - Location: line 38.
    - Evidence: serial runs exist only for (10,10)–(40,80) (FS and MONO as run) and for (20,20) and (40,40) (corrected MONO, FS with exact LU on S). The new meshes have none.
    - Fix: "for every solver and mesh run in serial".

17. **[LOW] The comparisons with the exact shifted inverse are inconsistent.**
    - Section 0 says "5-8 against 10-12" (line 35); Section 2.3 says "5-7 against 10-11" (line 224). Both mix left and right preconditioning.
    - Like-for-like at (20,20)/(40,40):
      - left: 6 and 7 against 12 and 12 (F_fsdirect);
      - right: 6 and 8 against 11 and 10 (Set B).

18. **[LOW] Misleading reason given for why the command-line fix works.**
    - Location: line 236, "(the key is not in the solver dictionary, so the command line works)".
    - Evidence: command-line options override `solver_parameters` regardless (petsctools/options.py:454-462).

19. **[LOW] Section 2.5(h) is stale, and results.csv with it.**
    - H_mono_monitor finished at 09:47: 100/100/100, rc 0, 635 s. results.csv still records it as "100".
    - H_mono_mumps_np16_nx80_nz80 is still running at 10:13: 17, then 100 (cap). This matches the PETSc-LU run.
    - Fix: re-run parse_table5.py and record_tables.py once it finishes.

20. **[LOW] Two side-finding numbers in Section 2.9.**
    - Line 478 "6e-10 for the old MONO": it is 6.2e-10 at 16×4 and 9.2e-9 at 32×8.
    - Line 471 "5.5e-6": no saved file contains it (analyse_nx4.txt is rotation-off, and dump_fs_rot_nx4.npz is not on disk). It cannot be verified.

21. **[LOW] Table 4 sentence and one Table 4 value (the sentence is touched by the diff).**
    - ex_article.tex:3737-3738 says "wall-clock cost of a single preconditioner application". The values are cumulative over the 14 PCApply calls of the timed solve (data_firedrake26_reordering/log_view_*: count 14; e.g. 1.6845 s / 14 = 0.12 s per application).
    - The (20,40) no-RCM entry 35.01 s should be 35.00 s (log: 3.5000e+01).
    - The R² sentence still says "Time(PCApply)" after the header rename. R² = 0.998 itself checks out.

22. **[VERY LOW] Three cosmetic points.**
    - fmg.c:6-32 shows the FMG cycle, not that restriction is the transpose of interpolation. That comes from mg.c:1108 → mgfunc.c:362, plus mat/interface/matrix.c:9199.
    - run.log prints "C1 0.1" (LB_time_slice.py:92 uses δ·dt^1.5), while the record and paper use C_1 = δ·dt = 1e-2 (ex_article.tex:2996). A one-line note would avoid confusion.
    - The paper text says "2.1×10^6" unknowns while the table says "2.05M" (2050560).

**CHECKS PASSED**

- **Set T (new Table 5):**
  - All 45 runs: counts, medians and [min, max] of the three repetitions, s/it and KSPSolve match the record. The paper tabular is byte-identical to table5_tabular.tex.
  - Counts are identical across repetitions.
  - Scatter: maximum 9.7 % (FS at (40,40)), median 3.7 %.
  - Per-application cost ratio: 10.7–16.6 %, which rounds to 11–17 %.
  - The (80,80) cumulative difference is 0.95 %; the (80,160)→(160,160) scaling factor is 2.11; dx = 469 m, dz = 25 m.
- **Unknowns and patch sizes:**
  - Unknowns 32320 / 128640 / 513280 / 1025280 / 2050560 equal the ksp_view row counts.
  - The formula 20nN_z + 4n follows from utils.py: RT_2 n(8N_z+2), DQ1 u_y 4nN_z, W_θ n(4N_z+2), DQ1 p 4nN_z. It matches the audit's field sizes 544/256/288/256 at 16×4.
  - Patch sizes 30N_z+8 and 38N_z+8 match pc_star_view_patch_sizes in every run (608 … 6088, old 3048).
- **Other sets and the generated tables:**
  - results.csv matches the raw logs for all completed runs; only the then-running monitor run differs.
  - Section 2.5 equals results_record.md; Section 2.7 equals results_final.md.
  - Every count in Sets A, A′, B, C, F, H and Y matches the raw logs.
- **Old 2026-06-01 logs:**
  - 2 processes; L_x = 300 km, L_z = 4 km (ε = 1.3e-2).
  - GMRES left, restart 30, rtol 1e-8, atol 0; Chebyshev(6), ASMStarPC, RCM, PCMG full.
  - Counts (PCApply − 1): 12/13/13/13 and 5/6/7/7. PCApply times are as printed. The (40,80) monolithic patch is 3048.
  - Table 4 logs: 2 processes, same geometry.
  - The ratio 1.9–2.4 and the today/June ratio 1.04–2.09 check out.
- **Serial vs parallel:** iteration counts are identical in all three solves for every serial/2-rank/16-rank pair that was run.
- **Source references:**
  - Correct: mg.c:1999 (inside PCCreate_MG, which starts at line 1985), mg.c:935-943, mg.c:559-561, assembled.py:72, asm.py:58, asm.py:62-64, fieldsplit.c:771/906/2133-2136, LB_time_slice.py:106-107, :317 and :371-381.
  - The ksp_views agree: the old monolithic run's finest smoother uses J as Amat and A_δ as Pmat; monofix uses A_δ for both.
- **Agent figures traced to their logs:** q(2) = −0.0359, q(2)² = 1.29e-3; 18 vs 159 eigenvalues; 0.0165; 92–97 %; 2.4e-2; 1.3e-16; 2.5e-15; 1.0e-16; 6e-16 and 6e-13; 9e-7 to 3e-6.
- **Error messages:** the GMRES left + unpreconditioned-norm error and the (10,10) 16-rank MPI error are as described.
- **Paper build (scratch copy):** 54 pages, no undefined references, no overfull box in the edited region. HEAD's Table 4 overflowed by 71.28 pt.
- **Other paper numbers:** patches 20.8–21.0 % smaller; "6 and 8" and "11 and 10" correct.
- **Side findings:** `bool(Constant(0.0))` is True (ufl/core/expr.py:310-312), and the U_mean advection is added twice (utils.py:120-134).

Files are in /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_numbers/:
- myparse.py and myparse.json (raw-log parse)
- check_T.py (Table-5 statistics)
- cmp_csv.py (results.csv check)
- paper_build/ and paper_build_head/ (paper builds)
