# completeness agent report (workflow wf_68296924-398, 2026-10-19)

1. **"COMMIT2" placeholder; the Table-5 paper edits are not committed.**
   - Location: record line 47 ("Commit COMMIT2 (pushed)"). Also the Section 2.8 heading and the Section 3 Files row ("the Table-5 commit: Section 0").
   - Evidence: at 10:08 the paper repo shows `M ex_article.tex`, and HEAD = origin/main = 0778f01. Task 2 never asked for a push.
   - Fix, if the edits are committed and pushed: put the hash in all three places ("Commit <hash>, pushed to origin/main").
   - Fix, if they are not: write "Not committed (your decision); `git -C ../Patch-Preconditioner-for-Boussinesq diff` shows the edits" in all three places.

2. **Wrong percentage, in the paper and the record.**
   - Location: ex_article.tex:3811-3812 ("cumulative PCApply time of the fieldsplit Schur approach is lower by $11$--$17\%$"). Record lines 44-45 ("11-17 % cheaper"), 425, 431 and 445.
   - Evidence: 1 - fs/monofix = 1 - 1.75/2.05 = 14.6 %, 1 - 5.76/6.66 = 13.5 %, 1 - 44.60/49.36 = 9.6 %, 1 - 94.11/106.10 = 11.3 %. The 11-17 % figure is monofix/fs - 1 (17.1, 15.6, 10.7, 12.7 %).
   - Fix (paper): "is lower by $10$--$15\%$".
   - Fix (record line 44-45): "the monolithic approach needs 11-17 % more PCApply time (the fieldsplit approach is 10-15 % cheaper), except 1 % at (80,80)". Change lines 425, 431 and 445 the same way. Keep "one monolithic application is 11-17 % more expensive" (that one is correct).

3. **Paper Table 5 caption says the smoother runs on "all three levels".**
   - Location: ex_article.tex:3832-3834; record line 397.
   - Evidence: in T_fs_rep1_np16_nx20_nz20/slice3D.txt the coarse level 0 is "preonly + lu". Chebyshev with pc_star appears only on mg_levels_1 and mg_levels_2. The monofix runs are the same.
   - Fix: "(Chebyshev(6) smoothing with the column-star patch and RCM-reordered patch factorisations on the two finer levels of the three-level hierarchy, direct LU on the coarsest)". Same wording in the record.

4. **The MUMPS run of Set H is still running, but the record shows "-".**
   - Location: line 359, the Section 4 row at line 515 ("MUMPS-patch run: see Section 2.5 (h)"), and the H register row.
   - Evidence: runs/H_mono_mumps_np16_nx80_nz80 has no status or done file. It was started at 09:47:26 by run_table5_fixed.sh (pid 3503720). Its first solve logged "CONVERGED_RTOL iterations 17" at 09:51, and it was still running at 10:08.
   - Fix: wait, rerun record_tables.py and fill it in. Otherwise write "running at the time of writing (first solve 17 iterations, as with PETSc LU) [open]" in 2.5(h) and Section 4, and add it to the Open list.

5. **Wrong explanation of why serial and parallel counts agree.**
   - Location: lines 366-367 ("pinned at exactly 2 ... hash of the global index"). Related: line 198 ("exactly 2 ... [0.2, 2.2] on every level") and the register's multigrid row ("= [0.2, 2.2] here").
   - Evidence from slice3D.txt ("eigenvalue targets used"): for the as-run MONO the finest-level estimate depends on the partition.
     - (20,20): max 2.20854 serial, 2.20695 on 2 ranks, 2.20452 on 16 ranks.
     - (40,40): 2.20641 / 2.20218 / 2.20657.
     - (40,80): 2.20661 / 2.20773 / 2.21075.
     - The same holds for all C, H, Ap_monocheb and Y_mono runs.
     - At (160,160) the finest level is 2.20004 for fs/fspmat and 2.20001 for monofix.
   - Fix: "For FS and the corrected MONO the largest estimated eigenvalue is 2 on every level (2.00003/2.00001 at (160,160)), so their Chebyshev interval does not depend on the partition. For the as-run MONO the finest-level estimate does depend on it (e.g. (20,20): 2.0078 serial, 2.0063 on 2, 2.0041 on 16 ranks), and its counts are identical nonetheless." Change the register's "= [0.2, 2.2] here" to match.

6. **An audit-agent shift sweep was run but never reported, and it qualifies the "to within one" claim.**
   - Location: line 512 ("(to within one) the same counts"). Also the register's common shift row (delta = 1e-4 "for every run") and the agents row (line 275).
   - Evidence: agents/theory_blind/runs/E6_*/log.txt, (8,8) i.e. a 32 x 8 fine mesh, GMRES left. Counts for FS, monofix (`pc_use_amat 0`) and as-run MONO:

     | delta | FS | monofix | as-run MONO |
     |---|---|---|---|
     | 1e-3 | 17/22/22 | 15/24/24 | 6/6/6 |
     | 1e-5 | 6/6/6 | 6/6/6 | 3/3/3 |
     | 1e-6 | 5/5/5 | 4/4/4 | 3/3/3 |

     E7 also ran meshes (4,8) and (8,4).
   - Fix: add these to the register and to 2.5/2.6. Qualify Section 4: "to within one at delta = 1e-4 on every mesh; at delta = 1e-3 on the 32 x 8 mesh with GMRES left they differ by 2 (22 vs 24); not checked with right preconditioning".

7. **Wrong attribution of the code change.**
   - Location: lines 500-503 ("not changed by this study (the diagnostic flags ... were added earlier the same day, before your commit a9a20ee").
   - Evidence that this study added them: LB_time_slice.py mtime is 06:08:50, the same second as the first runner start (run_table5_attempt1.log) and after the 05:35 request. The `--mono_smoother` help text says "Diagnostic (Table 5 of the paper)". `git show a9a20ee -- LB_time_slice.py` shows these four flags plus the `--patch` flag of 2026-10-16.
   - The check itself passed:
     - `git diff LB_time_slice.py` is empty.
     - With no flags given the parameters equal d81b1b3's: the `mono_smoother` default 'gmres' is the old literal, `inner_rtol`/`pc_side`/`outer_ksp` only apply when given, and `--patch star` reproduces the old inline block including the RCM keys.
     - utils.py and Nonlinear_Boussinesq_Irk_SC.py are untouched (mtime 10-16).
   - Fix: "This study added four default-preserving diagnostic flags (`--mono_smoother`, `--inner_rtol`, `--pc_side`, `--outer_ksp`) at 06:08, before the first run. Your commit a9a20ee contains them together with the `--patch` flag of 2026-10-16. Nothing has changed since (`git diff` empty)."

8. **The experiment register misses runs and settings (Molin's explicit requirement).** I compared every run folder's cmd.sh with the register.
   - A' row: FS "right" is GMRES right (`--richardson` keeps GMRES); MONO "right" is FGMRES right (the monolithic default). Neither is stated. The four Ap_monocheb_left_rich runs (GMRES left plus `-firedrake_2_ksp_norm_type unpreconditioned`, all failed) are not in the row.
   - Y row: Y_mono_natural_right_np2_nx{20,40} (natural ordering, FGMRES right: 7/6/6 and 8/8/8) are missing from the row and from table (d). Y_mono_rcm_mumps_np2_nx80_nz80, stopped by hand per run_table5_pivot.log, is also missing.
   - H row: H_mono_monitor_right_np16_nx80_nz80 (as-run MONO, FGMRES right, with `-firedrake_2_ksp_monitor`) is not listed.
   - Agents row is missing:
     - the delta values 1e-3, 1e-5 and 1e-6;
     - the meshes (4,8) and (8,4);
     - the D_* and E3 runs with `--tmax 100` and a fixed Chebyshev interval of 0.2, 2.2;
     - the code_audit flag_* runs, made through a modified copy (lb_audit_copy.py) with Richardson(0.5), GMRES(1), LU-inner and fixed-interval variants;
     - the variants and left/right sides in numerics/fs_vs_mono.py.
   - Everything in sets A, B, C, F and T matches the register exactly.

9. **Wrong cross-reference and overstated stall.**
   - Line 230: "(Section 2.5 (f))" should be "(e)".
   - Line 515: "3 x 100 iterations with left (2 and 16 ranks)" should be "17/100/100 with GMRES left (2 and 16 ranks), 100/100/100 with FGMRES right (16 ranks)".
   - Line 35: "does not converge at all at (80,80)" should be "does not converge within 100 iterations at (80,80) (17/100/100 with GMRES left, 100/100/100 with FGMRES right)".

10. **The residual-monitor run was made but its result is not reported.**
    - Evidence (H_mono_monitor_right run.log), solves 2-3: residual 42.02, then 3.3e-4 after 5-6 iterations, then 2.0e-5 at 100 (relative 4.8e-7). Solve 1: 14.47, then 1.4e-6 after 16 iterations, then 2.9e-7 at 100 (relative 2.0e-8).
    - Fix: add to 2.5(h)/2.6(h): "fast initial reduction, then stagnation at a relative residual of about 5e-7, above rtol 1e-8".

11. **The paper is now internally inconsistent about the smoother.**
    - Location: record 2.8 item 6 says the Appendix "still holds".
    - Evidence: the Appendix (ex_article.tex:4354-4359) describes `ksp_type:richardson`, six sweeps with theta, and links that to "the comparison of Section 5.6 ... untuned on both sides". The new Table 5 caption says Chebyshev(6).
    - Fix: say so in 2.8 item 6 and add it to the decisions list.

12. **The shift choice is not presented as a decision.**
    - Evidence: Table 5 keeps delta = 1e-4 (C_1 = 1e-2). ex_article.tex:2996-2997 says delta = C_1/dt "throughout, with nominal value C_1 = 10^{-3}", i.e. delta = 1e-5 at dt = 100 s.
    - Fix: state in 2.7 that delta = 1e-4 was kept for comparability with the old table and is 10x the nominal shift. List a re-run at delta = 1e-5 as Molin's decision. The E6 sweep (item 6) gives FS 6 = monofix 6 at delta = 1e-5.

13. **Section 1.4 (line 125) points to replacement text that is not on disk.**
    - It says the "exact replacement text for each is in edits.json and the referee reports".
    - No report files exist in agents/{section_referee,figure_referee}.
    - edits.json has no text for the "Section 5.1 / figs. 4-7 FGMRES, right, ten orders" item or for the Appendix-Richardson item.
    - It mixes entries already applied (A1, B10, B11, part of B8) with unapplied ones.
    - Fix: list the unapplied IDs (B1-B7, the rest of B8, B9, C1, C2) and say that no replacement text exists for the two items above.

14. **Section 0 serial claim is overstated** (line 38, "for every solver and mesh").
    - Serial runs exist only for (10,10)-(40,80) (Set A, FS and as-run MONO) and for (20,20)/(40,40) (F monofix and fsdirect). (10,10) cannot run on 16 ranks, and there is no serial run at (80,80) or above.
    - Fix: "for every solver and mesh run in serial".

15. **Section 0 "Other errors" (lines 39-41) omits two count errors that Section 2.1 shows.** FS at (40,40) and (40,80) was printed 12 but the logs say 13. MONO at (20,20) was printed 7 but the logs say 6.

16. **Machine description is wrong** (lines 260 and 516, "32-core").
    - lscpu: 2 x Xeon Silver 4215, 8 cores each, 2 threads per core, i.e. 16 physical cores and 32 hardware threads, 62 GiB. 16 ranks fill every physical core.
    - Fix: "16-core (32-thread) workstation".

17. **Missing label.** Sections 2.5 and 2.6 jump from (f) to (h). Rename (h) to (g) there, in Section 4 and in record_tables.py.

18. **An edit to the previous record is not listed.** claude_small_patch_robustness.md (mtime 06:49) lost its Section 8 (about 37 lines, the morning draft that is in commit 0b68318), replaced by a pointer. Add it to the Section 3 Files table.

19. **Loose wording on timings.**
    - Line 326 calls the Set Y timings "not comparable" but then states "4-6 x slower". The measured PCApply ratios are 2.7-6.8x (FS 2.7x/6.1x, MONO 3.9x/6.8x), taken under load.
    - Line 335 "takes 40-100 min" is an estimate. The measured (80,80) runs took 8 min on 16 ranks and 20 min on 2.

20. **"fspmat timing within the run-to-run scatter" (lines 399-400) is not true at the larger meshes.** At (80,80) fspmat is +2.1 % against a scatter of 0.3-0.4 %; at (80,160) it is +2.5 %. Say "within 2.5 % of fs".

21. **Unsupported number.** The rotation figure "relative 5.5e-6" (line 471) appears in no file under agents/code_audit (only the dump log, no analysis output). Save the output or drop the number.

22. **Confusing log label.** Every monolithic run.log prints "monolithic MG ASMExtrudedStar (star line patch)", but the ksp_view shows firedrake.ASMStarPC. Add a one-line note to the register.

**Checks that passed:**
- **Section 1 against `git show 0778f01`:** every item is accurate.
  - The commit holds the tex and 12 PNGs; figs. 4-12 point to `Image_restyled/`.
  - Section 5.7 is two paragraphs plus figs. 14 and 15; the Section 3.4 and Discussion edits, the caption edits for figs. 6, 10, 11 and 12, and the Section 5.1 missing-word fix are all there.
  - The script changes at 07:21:21 come before the PNG regeneration at 07:21:23-29 and the commit at 07:24. origin/main = 0778f01.
  - Scratch builds: 53 pages, no undefined references, and no new overfull boxes compared with eeda07e (two were removed).
- **Working-tree paper build:** 54 pages, no undefined references. The 71 pt (Table 4) and 25.7 pt (old Table 5) overfull boxes are gone, with none added.
- **2.8 items 1-5** match the uncommitted diff.
- **Tables:** the Section 2.5 tables equal results_record.md. Table 2.7 equals results_final.md, table5_tabular.tex and the paper's Table 5. The cost (11-17 % more per application), scatter (max 9.7 %, median 3.7 %), unknowns and dx/dz figures are correct.
- **Common parameters:** they match the code defaults (refinement 2, degree 2, rotation off, N = 0.01, initial buoyancy, boundary conditions, PCMG full/V, coarse LU, Chebyshev(6) under `--richardson`, star patch construct_dim 0, restart 30, null space).
- **Cited source lines:** mg.c:1999 and 935-943, assembled.py:72 and asm.py:58/62-64 say what the record claims. The as-run finest-level smoother uses Mat (firedrake_2_), i.e. J. MUMPS really took effect in Set Y.
- **June logs:** 2 processes, 300 km x 4 km, same host, `--rtol 1e-8 --maxit 100 --timing`.
- **Molin's questions:** the serial test is reported (Set A serial for both solvers, F monofix serial). The doubt is answered clearly (mechanism, fix, equal counts up to 2.05M unknowns).

Scratch files (two diffs, three build logs and .aux files) are in /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_completeness/.
