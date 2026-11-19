# Boussinesq_implicit

Firedrake solvers for the linear and nonlinear Boussinesq equations in a vertical slice. They are
implicit in time and use the shifted auxiliary-operator (Schur-complement) preconditioner with
multigrid and patch smoothers. This is the code behind the paper *Patch-Preconditioner-for-Boussinesq*.
Running it needs Firedrake with Irksome, for example the venv `~/developer-firedrake-2026-07/venv-firedrake-07`.

The repository was reorganised on 2026-11-15. Files were only moved, and paths were updated
where needed; nothing was deleted. The old → new path table and the full record are in
`../README.md` and `../tidy_up_2026-11-15/`.

## Layout

| Path | What |
|---|---|
| `utils.py` | Shared code: UFL residual forms (LB, SLB, nonlinear, Irksome variants), function spaces, the 3D-embedded extruded mesh hierarchy, constants (N² = 1e-4) |
| `LB_time_slice.py` | **Linear solver** used for the paper: implicit midpoint, shifted Schur-complement aux PC, MG with star/vanka patches; also `--monolithic`, `--U_mean`, timing and diagnostic flags |
| `Nonlinear_Boussinesq_Irk_SC.py` | **Nonlinear solver** used for the paper: Irksome GaussLegendre(1) + Newton, `IRKAuxiliaryOperatorPC` with the shifted Schur complement |
| `Nonlinear_Boussinesq_Mountain.py` | Test case: non-hydrostatic / `--hydrostatic` flow over a mountain (notes: `docs/Nonlinear_Boussinesq_Mountain_summary.md`) |
| `Nonlinear_Boussinesq_Lock_Exchange.py` | Test case: lock-exchange gravity current, experimental (notes: `docs/Nonlinear_Boussinesq_Lock_Exchange_summary.md`) |
| `claude_paper_shift.py` | The analysis-derived shift δ = K δ_thr(dt, ε, N), used by the K-sweep generator, the plotter and the study scripts |
| `plot_ar.py`, `plot_dx.py`, `plot_dz.py`, `plot_shift.py` | Plots for the paper's sweeps (δ = C1/dt); run inside a data folder |
| `plot_claude_paper.py` | Plots and `summary_iterations.csv` for the K sweeps; run inside a data folder |
| `run_claude_paper.sh`, `claude_paper_run_one.sh` | Parallel sweep driver and its per-job runner |
| `Bash_Generator/` | Command generators that print the sweeps: `bash.py` / `bash_nonlinear.py` (the paper's original linear / nonlinear sweeps), `bash_paper_patch.py` (the same sweeps with `--patch star\|vanka`), `bash_claude_paper.py` (K sweeps), `Colin_bash.py` (unrelated template) |
| `sweeps/<category>/test_*.sh` | Serial sweep scripts, each a list of `mpiexec` commands plus a plot command, sorted by the data folder they produce: `paper_robustness/` (26), `claude_paper/` (30), `early_prototypes/` (1) |
| `docs/` | Notes: `Bug_advection_operator.md`, the mountain and lock-exchange summaries; `docs/claude_studies/` has the four Claude study records (K-shift tests, column patch, small patch / pcfix, paper update and Table 5) |
| `claude_column_patch/`, `claude_small_patch/`, `claude_patch_gmres6/` | Claude studies with their code and results (star vs single-column patch; small patch, pcfix diagnostics, Table 5 timing scripts, paper figure scripts; GMRES(6) robustness study). Their run folders, job lists, logs and agent scratch are in `../data/{column_patch,small_patch,patch_gmres6}/` since 2026-11-18; on the workstation the old paths are symlinks to them (ignored by git) |
| `archive/` | Superseded code, kept unchanged: `old_solvers/` (LB_Irk, LB_time_slice_separate_save, Monolithic_Boussinesq, Nonlinear_Boussinesq_slice, Nonlinear_Irk_Boussinesq_direct, SLB), `small_experiments/` (mesh-checkpoint MFEs, RCM reordering), `Colin_toy/`, `Flip_MFE_test/`, `Flip_SLB_test/`, `LB_test/`, `Plot/` |

## Where the data goes

The scripts write their output into the **current folder**. By convention, every run happens
inside a data folder `yr2/data/<category>/<name>/` next to this repository (categories:
`early_prototypes`, `nonlinear_development`, `paper_robustness`, `meanflow`,
`firedrake26_checks`, `mountain`, `claude_paper`, `patch_gmres6`; the bookkeeping of the Claude
studies is in `small_patch` and `column_patch`). Seen from such a folder, the code is at
`../../../Boussinesq_implicit/`.

Do not start a solver, test case or serial sweep from inside this repository; its output would
land here. The batch runners `sweeps/paper_robustness/test_script.sh` and
`test_nonlinear_script.sh` are the exception: they `cd` into their data folders themselves.

## Running

```bash
# one run: make (or enter) a data folder, then call the solver from there
mkdir -p ../data/<category>/<new_run> && cd ../data/<category>/<new_run>
mpiexec -n 2 python ../../../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 \
    --dt 50.0 --tmax 100.0 --shift 2e-05 --maxit 100 --dt_test --richardson --reordering

# a whole serial sweep, from inside its data folder (the script ends with the plot command)
bash ../../../Boussinesq_implicit/sweeps/<category>/test_<...>.sh

# parallel sweeps, from this folder:
#   data -> ../data/<category>/<prefix>_<mode>_<test><suffix>/, serial form -> sweeps/<category>/,
#   finished jobs (same command) are skipped
./run_claude_paper.sh -m <linear|nonlinear> -t "<dt ar dx dz>" [-s <suffix>] [-x "<generator args>"]   # K sweeps -> data/claude_paper/
./run_claude_paper.sh -g Bash_Generator/bash_paper_patch.py -D data -S test_paper ...                              # paper sweeps -> data/paper_robustness/
#   -n: only write the command lists and job folders; -p: only re-plot; -C <category>: another data/ sub-folder
#   (the exact commands used for the existing sweeps are in claude_small_patch/run_all_sweeps.sh and
#    docs/claude_studies/claude_paper_tests_README.md)

# paper figures: they read data/paper_robustness/ and data/claude_paper/ and write into
# ../Patch-Preconditioner-for-Boussinesq/Image_restyled/ or Image_patch/ (overwriting the figures there)
python claude_small_patch/restyle_paper_figures.py --set paper
python claude_small_patch/paper_fig_single_column.py
python claude_small_patch/restyle_fig14.py
python claude_column_patch/paper_figure.py
```

Some archived scripts `import utils` from the repository root. To run one of those, from this
folder: `PYTHONPATH=$PWD python archive/old_solvers/SLB.py`. Several no longer run with today's
`utils.py`; see `../README.md`, section 5.

## Notes

- The data folders made by the driver contain `commands.txt` and `jobs/*/cmd.sh`. These record
  the commands as they were run, with the old `../Boussinesq_implicit/` (or absolute) paths.
  The driver regenerates `commands.txt` on every call.
- Re-running the original serial sweeps (`test_{dt,ar,dx,dz}.sh`, `test_nonlinear_*.sh`) at today's
  HEAD does not exactly reproduce the paper's data, because some solver defaults have changed since.
  The `test_paper_*_star.sh` scripts are the like-for-like re-runs.
- `*.png`, `*.txt`, `*.out`, `*.log`, `*.vtu`, `*.pvd`, `*.h5`, `*.npz` and files named `status`, `done` or
  `finished` are git-ignored everywhere in the repository, as are the symlinked run folders of change-log item 5.
  A new `.npz` or `.png` result therefore needs `git add -f`.

## Paper repository updates (2026-11-18)

Changes made in `../Patch-Preconditioner-for-Boussinesq/` (uncommitted; `git diff` there shows them). Only
`ex_article.tex` and `ex_supplement.tex` changed; no figure file was touched.

**1. Check of the paper's numbers against the workstation reruns** (`claude_patch_gmres6/`, same solver configuration
as the paper: GMRES(6) level smoother, FGMRES right, back-substitution, 2 ranks).
- Linear (paper §5.2, Table 4, Fig. 2; supplement comparison): all reproduced.
  - Mesh sweeps at Δt = 50 s: 6 or 7 iterations.
  - Aspect-ratio table: 17, 12, 7, 4, 3, 2, 2.
  - Shift-sweep failure set and the quoted counts (46, 54, 29; 2–13 at Δt = 10⁴ s).
  - Single column within one iteration in 152 of 185 runs; 33 runs more; the same 14 failures; 11 against 2 at
    Lz = 250 m, Δt = 10³ s.
  - Differences of ±1 iteration only, all at noisy points (parallel MUMPS coarse LU): C1 = 1, Δt = 100 s gives 85
    here against the paper's 86; at C1 ≤ 1e-4 the single column needs at most 7 against 2 here, against the paper's 6.
  - Two statements that hold for one run but not another were reworded so that they hold in both data sets (item 2).
- Nonlinear "6 iterations with the linear equations at C1 = 1e-4, Δt = 50 s" (§5.3): reproduced (6, both patches,
  `--linear_forms`; `data/patch_gmres6/paperconv/linear_forms_check/`).
- Nonlinear sweeps (§5.3, §5.4): checked on 2026-11-18 against 440 fresh runs in the paper's shift convention
  (`--shift C1/(2Δt)`; both patches; all exit 0). Scripts: `claude_patch_gmres6/make_paperconv_jobs.py`,
  `claude_patch_gmres6/analyse_paperconv.py`. Data: `data/patch_gmres6/paperconv/`. Results:
  `claude_patch_gmres6/figures/80_paperconv_nonlinear/{counts.csv,checks.txt}`.
  - Confirmed:
    - 8 or 9 iterations in 2 Newton steps at Δt = 50 s;
    - aspect ratio 7, 7, 6, 5, 4 at Δt = 100 s;
    - C1 ≤ 1e-5 lowest at every Δt; C1 = 1e-4 gives 21–35; level-off at 3–7; 21 against 8;
    - Newton 2–3 steps; the same 2 failures for both patches;
    - 106 admissible runs with median ratio 2.7; single-column Nz medians 8 → 12 against the star's 3 → 4;
      112 other runs with median ratio 1.0.
  - Differences, reworded so that they hold for both data sets:
    - "25–120" (C1 ≥ 5e-4): 121 here → "25 to about 120";
    - "at most 6.3" (admissible ratio): 5.7 here → "at most about 6";
    - "equal in 63": 62 here → "equal in more than half".
  - Repeatability: these runs and the earlier doubled-shift runs share 29 identical commands per patch, and give
    the same counts in all 29, for both patches.
  - The repeatability sentence of §5.4 now reads: through the parallel coarse direct solve, repeated runs can differ
    by up to three iterations in single-column counts at large Δt and by one near the cap; all other counts are
    reproducible. Basis: linear 11–14 over ten runs; star 85/86 and 120/121 between the Mac and the workstation.

**2. Single column as an equal option in the main text** (decision of 2026-11-18: option 3, the article may exceed
28 pages; it is now 30, the supplement 11).
- New article §5.4 `sec:patch-comparison` "Star and single-column patches", placed before "Flow over a ridge".
  - It takes over the experimental part of supplement SM2: patch sizes, the linear and nonlinear comparison, the
    cost sentence, and the two figures `fig:single-column-linear` / `fig:single-column-nonlinear`
    (`Image_restyled/single_column_{linear,nonlinear}.png`, unchanged).
  - The SM2 TODO about the 30 % moved with it.
  - Rewordings against SM2:
    - "fewer only once (85 against 86)" → "where the single column needs fewer, it is by one iteration";
    - "all but one (94 against 91 at Nx = 10, Δt = 1 s)" → "apart from runs at Δt = 1 s in which both patches need
      more than 75 iterations";
    - "at most 6 against 2" → "at most 7 against 2".
  - New sentence: star counts are reproducible; single-column counts at large Δt vary by up to three iterations
    between repeats, through the parallel coarse direct solve.
  - New closing paragraph: both patches are robust in resolution and aspect ratio; the star needs fewer iterations at
    small shifts, the single column has smaller patches (about 30 % cheaper per iteration, still the TODO number).
- §3.4 `sec:inner-solvers`:
  - "We therefore use the star patch throughout" is replaced by: the single column is an alternative for k ≥ 2;
    §5.4 compares the two patches (it converges wherever the star does, with smaller patches but more iterations at
    small shifts, the more so the flatter the cells); the other experiments use the star.
  - The summary sentence now reads "one multigrid cycle with star or single-column patches".
- §5.1 `sec:solver-config`: the paragraph with the single column's 30 %, 5.5× and 2–3× was replaced by the solver
  setting of the single column (`firedrake.ASMVankaPC`, `pc_vanka_construct_dim:1`, `pc_vanka_column:1`, all
  other options unchanged) and a pointer to §5.4.
- Introduction: one sentence that in a slice with k ≥ 2 single columns plus the coarse grid also capture the kernel,
  and that the experiments compare both patches. The outline now names the supplement's "kernel argument for the
  single-column patch" instead of "a comparison of patch sizes".
- Conclusions: "The star patch is used throughout; … a cheaper alternative where its iteration penalty is small" is
  replaced by:
  - the single column is an equally robust alternative in the slice for k ≥ 2;
  - it costs 2–3× the star's small counts at admissible shifts in the nonlinear problem;
  - the star is required for k = 1 and remains the choice in 3D.
- Supplement:
  - SM2 is renamed "The single-column patch" (label `sm:single-column` kept) and now holds only the definition, the
    patch sizes and the kernel argument (k ≥ 2 sufficient, k = 1 fails), with a pointer to §5.4.
  - `app:petsc-options` got one sentence on the single-column options (`ASMVankaPC`, `pc_vanka` prefix,
    `construct_dim 1`).

**Unchanged and still open**
- Solver settings: the article and supplement already matched the code (GMRES(6), FGMRES restart 30, full MG,
  RCM patch LU, direct coarse solve, back-substitution, EW 1e-2, atol 1e-9, maxit 100); nothing else was changed.
- ~~The red timing TODOs (§5.6 Cost, the 30 %)~~: done the same evening, see item 3.
- Build check (scratch copy, `latexmk`): article 30 pages, supplement 11; no undefined references; no new overfull
  boxes.

**3. Workstation timing run and the two red TODOs (2026-11-18, evening)**
- Run: `claude_small_patch/table5/run_timing_tables_gmres6.sh`, 64 runs one at a time, 18:56–20:49, machine otherwise
  idle (sar: 6.4 % CPU during the 2-rank set = our 2 ranks). Runs in `../data/small_patch/table5_runs_gmres6/`
  (symlinked as `claude_small_patch/table5/runs_gmres6`, item 5), log `claude_small_patch/table5/run_timing_gmres6.log`,
  summary `results_gmres6.{md,csv}` (made with the Mac's revised `parse_timing_gmres6.py`
  of origin/main 24b7f6a, run from a scratch copy; the local copy is the older one).
- Checks: 64 of 64 present, rc 0, no problems; iteration counts identical in all repetitions; T_fs (20,20) 11,
  B_fs (20,20) 11, set C 8 / 10 applications, all as on the Mac.
- PCApply spread over the three repetitions: 0.7–9.1 % for the 16-rank set T (October: at most 9.7 %), but 17–37 %
  for the 2-rank set R. Not other load (see sar above); most likely the unbound placement of 2 ranks
  (`--bind-to none`) on this 2-socket machine. The paper uses medians; the trends (linear with RCM, superlinear
  without) hold for the minima and maxima as well. If wanted: rerun set R with `--bind-to core` (about 40 min).
- Set C (nonlinear solver, (20,20), 2 ranks): per application star 0.3228 s, single column 0.2660 s (medians):
  17.6 % cheaper (the three paired repetitions give 14 %, 30 %, 18 %); break-even factor 1.21. Recount with
  `counts_gmres6.csv` (the paper's data): the single column is faster in 81 of 112 runs above the admissible range
  and 3 of 106 inside (October's 30 %: 97 and 6). The workstation reruns (`80_paperconv_nonlinear/counts.csv`)
  give 81 and 4.
- Paper edits (`ex_article.tex`; both TODOs removed; the supplement needed no change: its RCM and "untuned"
  sentences still hold):
  - `tab:rcm`: new rows (bandwidths unchanged); caption "the 10 or 11 preconditioner applications of one FGMRES
    solve … medians of three runs" (was "14 … of one GMRES solve"); text R² 0.998 → 0.992, and "RCM reduces
    PCApply fourteenfold at N_z = 80, from 128.3 s to 9.2 s" (was "reaching 139.7 s").
  - `tab:monolithic`: new rows; iterations now identical for both solvers (11, 10, 10, 9, 9); caption
    Chebyshev(6) → GMRES(6).
  - Text: 9 to 11 iterations; cost linear in the unknowns (5.6e-5 → 4.4e-5 s per unknown); exact solve 11 and 10,
    "exactly as many as with one cycle". Comparison paragraph: the Chebyshev-interval argument is gone; the two
    cycles need the same iterations at every resolution, and the fieldsplit preconditioner is 6–17 % cheaper (was
    "monolithic 11–17 % higher except (80,80)").
  - §3.3: "Section 5.6 shows that the two cycles need the same number of iterations …"; first sentence "On each
    patch the two routes relax the velocity in the same way".
  - §5.4: 30 % → 18 %, factor 1.4 → 1.2, 97/112 → 81/112, and "the star patch is faster in 103 of the 106 runs".
  - Conclusions: "about 10–15 % cheaper … at four of the five resolutions" → "6–17 % cheaper … at the five
    resolutions tested"; the single column's "faster in only 6 of 106" → the star "in 103 of 106".
  - Build (scratch copy): article 30 pages, supplement 11, no undefined references, same overfull boxes as HEAD.
- Solver code: local HEAD 6c4b0cf defined `--mg_smoother` twice in `LB_time_slice.py` (merge), which crashed every
  run. I deleted the duplicate line locally (uncommitted). origin/main 24b7f6a has the Mac's equivalent fix: for
  every `--mg_smoother gmres` run both versions give the same solver (outer FGMRES). Not pulled; to sync:
  `git checkout -- LB_time_slice.py`, then `git pull --ff-only`.

**4. Wording revised along the anti-defensive principles (2026-11-18, at Molin's request)**
Only passages I wrote; every number is exact for the paper's own (Mac) data. Two independent checks followed: one
for accuracy (no overclaim), one against the principles (`~/.claude/skills/anti-defensive-writing-en/SKILL.md`).
- §5.4: removed the repeatability sentence (parallel coarse solve) and the list of the single column's larger counts
  (5.5×, 7 against 2, 5 against 3). Now: "The other runs are at Δt ≥ 100 s, where the star patch needs two or three
  iterations and the single column four to eleven, and at Δt = 1 s, where both need more than 75"; "In the
  nonlinear problem both patches converge in the same 218 of the 220 runs"; exact numbers where item 1 had hedged:
  "25–120" (not "25 to about 120"), "equal in 63" (not "more than half"); "at most about 6" dropped. Inside the
  admissible range "the star patch profits from the small shift: three or four iterations at every N_z (medians),
  against eight at N_z = 10 and twelve at N_z = 80 for the single column". Closing: "Both patches converge in the
  same runs, and the choice between them is one of cost"; the single column "is the faster patch wherever it needs
  fewer than about 1.2 times the star's iterations, as in 81 of the 112 …"; "the star patch is faster in 103 of
  the 106 runs and is the patch of choice".
- §3.3: "With star patches, the two routes have equivalent patch problems" (the argument holds for star patches
  only); "Section 5.6 shows that the two cycles need the same number of iterations"; no caveat that the two GMRES
  cycles are not identical. Cost comparison: "The two cycles need the same number of iterations at every
  resolution, and they differ in cost: … the fieldsplit preconditioner is 6–17 % cheaper".
- §5.6: "its cost per layer falls from 0.185 to 0.115 s; without it, the cost per layer grows almost fivefold, and
  at N_z = 80 RCM reduces PCApply fourteenfold"; "its cost grows linearly with the number of unknowns"; "With A_δ
  solved to a relative residual of 10⁻¹⁰ in place of one multigrid cycle, FGMRES needs 11 and 10 iterations …,
  exactly as many as with one cycle"; tab:monolithic caption "six GMRES iterations preconditioned by the star
  patch"; tab:rcm (20,20) 3.454 s, as in `results_gmres6.md` (the median 3.4545 is a rounding tie).
- §3.4: "with smaller patches but more iterations at small shifts, the more so the flatter the cells" → "with
  smaller patch problems, and the star converges faster at small shifts".
- Introduction: "the experiments compare the two patches" → "in the experiments they converge wherever the star
  patches do".
- Conclusions: "6–17 % cheaper than a monolithic patch smoother at every resolution tested, up to 2.05×10⁶ unknowns,
  with the same iteration count"; the single column is "an alternative to the star: it converges wherever the star
  does, and at Δt = 50 s, the time step of the mesh sweeps, its counts are as independent of the resolution as the
  star's"; "The star patch is the faster choice at admissible shifts in the nonlinear problem, in 103 of 106 runs".
- The first rewrite had gone too far in four places, which the accuracy check caught and which are now narrowed as
  above: "Both patches are robust in the resolution and the aspect ratio" (the single column's counts at admissible
  shifts rise 8 → 12 with N_z and as the domain thins), "the star … captures the kernel without the help of the
  coarse grid" (§3.4 and SM3: the star also needs the coarse grid for the uniform flow), "the experiments show that
  they are robust as well" and "a robust alternative".
- Typo fixed: "Star-patch" → "star-patch" (§3.3).
- New file `../Patch-Preconditioner-for-Boussinesq/numerical_results_summary_2026-11-18.md` (3994 characters): the
  answers to Molin's five questions (robustness, influence of the patch, why the single column works, how to read the
  figures, what is missing), fact-checked against the data.
- Not changed (not my text), for Molin to decide; the check against the principles suggests: the last paragraph of
  the Conclusions ("shown numerically, not proved", "do not exercise the worst-case pressure mode", "Only the
  experiments are confined to the vertical slice"), which closes the paper on limitations; "agrees qualitatively"
  (Conclusions, abstract); "The agreement is qualitative … not exercised" (§5.2); "Periodicity only spares the slice
  its lateral walls" (§5.1); "For Δt ≲ 2×10³ s this lies above the worst-case admissible range" (§5.1); "while the
  wave train is still developing" (§5.5); and "the solver converges on a terrain-following mesh", where §5.5 shows one
  Newton step with two FGMRES iterations per time step.
- Found by the checks, pre-existing: the no-RCM bandwidths of tab:rcm (531, 1051, 1571, 2091) equal 52 N_z + 11 for
  N_z = 10, 20, 30, 40, not for the table's 10, 20, 40, 80 (at N_z = 10 the star patch has only 308 unknowns, so 531
  must be the full bandwidth 2 max|i−j| + 1). Probably the last two rows should read 2091 and 4171; they need
  recomputing from the patch matrices. Also: `LB_time_slice.py --no_rotation` is `store_false`, so passing it turns
  rotation ON (no paper run passes it, so f = 0 everywhere as intended).

**5. Repository clean-up (2026-11-18, at Molin's request)**
Why: my studies had kept their run bookkeeping inside this repository, and commits that added everything had put
it on GitHub: 858 tracked files, 133 MB of them under `claude_small_patch/table5/agents/` (three 43 MB `dense.npz`
dumps and a 2 MB PDF among them). Nothing was deleted; the record and a tested undo script are in
`../tidy_up_2026-11-18/` (section 7 of `../README.md`).
- Moved to `../data/` (16 moves, 3577 files, checked against before/after manifests and by two independent audits):
  - `claude_small_patch/table5/{runs,runs_gmres6,agents}` and `claude_small_patch/diag_item3/runs` →
    `../data/small_patch/{table5_runs,table5_runs_gmres6,table5_agents,diag_item3_runs}/`;
  - the run, scratch and log folders of `claude_column_patch/` (`timing/runs`, `code_verification/runs`,
    `firedrake_inner_solve/{scratch,logs}`, `firedrake_algebraic/work`, `theory_kernel/{solver_runs,work}`,
    `toy_model/logs`, `my_checks/degree1_runs`) → `../data/column_patch/<folder>_<subfolder>/`;
  - `claude_patch_gmres6/logs` → `../data/patch_gmres6/logs/`, and `claude_patch_gmres6/joblist_main.txt` into it;
  - `sweeps/patch_gmres6` (generated serial command lists) → `../data/patch_gmres6/sweeps/`.
- At each old path except `claude_patch_gmres6/joblist_main.txt` (now `claude_patch_gmres6/logs/joblist_main.txt`)
  a relative symlink to the new place, so every script's default path (`runs/`, `logs/`, `scratch/`, `work/`,
  `solver_runs/`, `runs_gmres6/`) still works on the workstation. The symlinks are ignored by git.
- What went along besides run records (`cmd.sh`, `status`, `done`, logs, outputs) and the 24 generated sweep lists:
  only scratch scripts: the agents' own 36 scripts, two exploration scripts of `claude_column_patch/firedrake_algebraic/work`
  and two `run.sh` of run folders. No kept code calls them; `docs/claude_studies/claude_paper_update_and_table5.md`
  cites four agent scripts, and its note points to the new place. (`code_verification/runs/parallel_patches/run.sh`
  reaches `02_patch_dofs.py` through `$HERE/../../`, which now points into `../data/column_patch/`; re-running it
  needs an absolute path.)
- Kept in the repository: the study code (including `claude_small_patch/diag_item3/Nonlinear_pcfix.py`, whose line
  ranges `claude_notes` executes), write-ups, results tables and small figures (the `.png` files on disk only:
  `*.png` is git-ignored). The two results tables of the column-patch timing study, `timing_summary.md` and
  `timing_all_runs.csv` (quoted in `claude_column_patch/paper_review_notes/numbers/numbers_check.md`), were copied
  from the moved run folder to `claude_column_patch/timing/`; git records them as renames.
- `.gitignore`: the 15 symlinked paths, plus `status`, `done`, `finished` (run-folder markers) and `*.npz`.
- Staged now (not committed): the 856 deletions, the 2 renamed tables, `.gitignore` and the notes in
  `docs/claude_studies/`. A plain `git commit` therefore takes the deletions together with the ignore rules. The
  path updates in this README and in `claude_patch_gmres6/{README,PROCEDURE}.md` are unstaged, together with these
  files' earlier changes.
- To finish: the sync is done (item 6); the staged state survived it. The commits and the push are left to Molin:
  the Claude Code auto-mode safety check did not let me commit the 856 deletions.
  The files leave GitHub with the push and stay in the history. When the Mac pulls that commit, git deletes the Mac's
  copies of the 858 files (they only ever came from GitHub; the originals are in the workstation's `yr2/data/`);
  anything the Mac made itself in these folders is untracked there and stays, and the pull stops if the Mac has
  local edits to any of the 858.
- `claude_small_patch/table5/WORKSTATION_TIMING_RUN.md` (the version arriving with the pull) suggests
  `git add -f claude_small_patch/table5/runs_gmres6 …` for a commit of the timing runs. Do not: that path is now an
  ignored symlink, and `-f` would commit the symlink. Commit `results_gmres6.{md,csv}`; the runs stay in `../data/`.
- Paths updated in this README, `claude_patch_gmres6/{README,PROCEDURE}.md`, and a note at the top of the three
  `docs/claude_studies/` records that refer to the moved folders (their old paths still resolve on the workstation).

**6. Table 2 re-run with pinned ranks, bandwidth correction, LSR, sync (2026-11-18/19, Molin's decisions)**
- Table 2 (set R) re-run with the two ranks pinned, one per socket (`mpiexec --bind-to core --map-by socket`;
  runner `../data/small_patch/table5_runs_gmres6_R_bound/run_timing_R_bound.sh`, a copy of
  `run_timing_tables_gmres6.sh` that changes only the binding; runs in that folder, summary
  `claude_small_patch/table5/results_gmres6_R_pinned.{md,csv}`). 24 of 24 runs rc 0, the same iteration counts
  as before (11, 11, 10, 10). The repetitions now agree to within 4% (11% for no RCM at N_z = 10), against 17–37%
  unpinned. The unpinned medians were too high at RCM (20,10) and (20,20) and no RCM (20,40), by up to a third
  (their minima match the pinned medians); the paper and the LSR now use the pinned medians: RCM 1.560, 2.870,
  4.874, 9.273 s; no RCM 2.749, 9.460, 32.08, 127.3 s; affine fit 0.56 s + 0.109 s per layer, R² = 0.999.
  An LTeX language server (VS Code) ran on the same socket during repetitions 2–3; those agree with the clean
  repetition 1.
- No-RCM bandwidths of Table 2 corrected to 531, 1051, 2091, 4171 (52 N_z + 11) in the paper, the LSR,
  `parse_timing_gmres6.py` and `WORKSTATION_TIMING_RUN.md` (which also no longer says to `git add -f` the runs).
- LSR (`Tex_Repo/LSR_Copied_from_ESA/main.tex`, uncommitted): all nine timing TODOs settled with the same numbers
  and wording as the paper, at the LSR's level of detail. Builds at 87 pages, no undefined references, the same
  overfull boxes as its HEAD. TinyTeX here lacks csquotes, algorithm2e, todonotes, babel-english and the Charter
  fonts; the build used them from a scratch tree (TinyTeX itself was not changed).
- Sync: the local `LB_time_slice.py` fix was discarded and `git pull --ff-only` brought 24b7f6a (the Mac's
  equivalent fix, the revised parser, `counts_gmres6.csv`).
