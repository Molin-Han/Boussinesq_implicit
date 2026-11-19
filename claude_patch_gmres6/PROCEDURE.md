# Procedure of the GMRES(6) star / single-column robustness test (2026-11-16)

This file covers how the test was run, what it shares with the earlier robustness tests, what differs, and where the
results are. The findings themselves are in `README.md` next to this file. All paths are relative to `yr2/`.
Since the clean-up of 2026-11-18 (`tidy_up_2026-11-18/`) the job lists and pool logs are in `data/patch_gmres6/logs/` and
the generated serial command lists in `data/patch_gmres6/sweeps/`. On the workstation
`Boussinesq_implicit/claude_patch_gmres6/logs` and `Boussinesq_implicit/sweeps/patch_gmres6` are symlinks to them, so
the commands below work unchanged.

## 1. Summary

- **What was run.** The paper's four robustness sweeps for the linear and the nonlinear code: time step (with the shift
  constant C1), aspect ratio, Nx and Nz, with the shift delta = C1/dt. Each sweep was run once with the star patch and
  once with the single-column patch. In BOTH codes the multigrid level smoother is GMRES(6) and the outer Krylov method
  is FGMRES with right preconditioning. These 854 runs used 2 MPI ranks each. On top of them:
  - 320 runs of 9 points per code on 1-16 ranks (the MPI study);
  - 10 runs that locate the source of the run-to-run noise;
  - 6 runs that diagnose the two runs that aborted on 16 ranks.
- **The same as before:**
  - the equations of the linear code;
  - the meshes and the exact parameter set of every sweep;
  - the shift rule;
  - the structure of the preconditioner;
  - 2 ranks per run;
  - iteration caps, tolerances and the Newton settings;
  - the output files.
- **Different:**
  - the smoother and the outer Krylov method: against the paper in both codes, against the 2026-10-16 rerun in the
    linear code only;
  - the patch: the paper used only the star;
  - the nonlinear advection term, corrected after the paper's runs;
  - the amplitude of the linear code's initial state, which does not affect its counts;
  - a few solver flags that are now passed explicitly;
  - the software stack;
  - the scheduling of the runs;
  - the added MPI runs, noise runs and analysis.

## 2. The tests that are compared

`<t>` = dt, ar, dx, dz; `<p>` = star or vanka (vanka = the single-column patch, `firedrake.ASMVankaPC`).

| | A. Paper sweeps | B. Star / single-column rerun | C. This test |
|---|---|---|---|
| runs made | linear 2026-05-19, nonlinear 2026-05-20 | 2026-10-16 | 2026-11-16 |
| used for | the paper's figs. 4-12 | supplement SM2 (patch comparison); `Boussinesq_implicit/docs/claude_studies/claude_small_patch_robustness.md` | `README.md` next to this file |
| data | `data/paper_robustness/data_linear_<t>_deg2/`, `data_nonlinear_<t>_deg2_ew/` | `data/paper_robustness/data_linear_<t>_deg2_<p>/`, `data_nonlinear_<t>_deg2_ew_<p>/` | `data/patch_gmres6/data_linear_<t>_deg2_gmres6_<p>/`, `data_nonlinear_<t>_deg2_ew_gmres6_<p>/` |
| command generator | `Boussinesq_implicit/Bash_Generator/bash.py`, `bash_nonlinear.py` | `Boussinesq_implicit/Bash_Generator/bash_paper_patch.py` | `bash_paper_patch.py`, plus `--extra=--mg_smoother=gmres` for the linear code |
| serial command lists | `Boussinesq_implicit/sweeps/paper_robustness/test_<t>.sh`, `test_nonlinear_<t>.sh` | `Boussinesq_implicit/sweeps/paper_robustness/test_paper_*_<p>.sh` | `data/patch_gmres6/sweeps/test_paper_*_gmres6_<p>.sh` |
| runs | 207 linear + 220 nonlinear (star) | 2 x (207 + 220) | 2 x (207 + 220) = 854, + 320 MPI + 16 diagnostic |

The nonlinear data folders of A also hold older files from 2026-05-02/03, from runs at shifts that are not in the
sweeps (e.g. `error_dt10000.0_shift1e-11.out`). The paper's figures use the runs of 2026-05-19/20.

## 3. What is the same in A, B and C

| item | setting |
|---|---|
| model | 2D vertical slice of the Boussinesq equations, f = 0 (no rotation), N = 0.01 1/s, no mean flow (U_mean = 0), degree-2 elements ("deg2"). Linear code `LB_time_slice.py` (one linear solve per time step); nonlinear code `Nonlinear_Boussinesq_Irk_SC.py` (Irksome, Newton). The initial buoyancy perturbation is the same bump at x = Lx/2 with half-width 5 km. The linear code solves the same equations in A, B and C. Later changes to `utils.py` only add optional terms (custom Coriolis parameter and N^2, mean-flow advection), which are zero or off here; for the amplitude of its initial state see Section 4. For the nonlinear code see Section 4. |
| mesh | periodic interval base mesh of Nx cells, Lx = 300 km; 2 uniform horizontal refinements (3-level hierarchy, 4 Nx columns on the finest level); Nz layers on every level (no vertical coarsening); height Lz. Nominal: Nx = Nz = 20, Lz = 4 km. |
| parameter sets | identical, checked by parsing the commands: every (Nx, Nz, Lx, Lz, dt, tmax, shift, maxit) tuple of A's test scripts appears in C and vice versa, for all 16 sweep/patch combinations. Linear: 207 runs per patch (dt sweep 54, ar 63, dx 45, dz 45). Nonlinear: 220 runs (70, 50, 50, 50). |
| parameter lists | **linear:** dt = 1, 5, 10, 50, 100, 500, 1000, 5000, 1e4 s; dt sweep C1 = 1, 0.1, 0.01, 1e-3, 1e-4, 1e-5; ar/dx/dz sweeps C1 = 1e-3; Lz = 12, 8, 4, 2, 1, 0.5, 0.25 km. **nonlinear:** dt as linear plus 5e4 s; dt sweep C1 = 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6; ar/dx/dz sweeps C1 = 1e-5; Lz = 12, 8, 4, 2, 1 km. **Both:** Nx, Nz = 10, 20, 40, 60, 80. |
| shift | delta = C1/dt, passed as `--shift` (the paper's design) |
| what a run does | `--tmax 2 dt` with the sweep's `--<t>_test` flag. The same solve is done twice: the first gives the reference solution, the second is monitored (error history and iteration counts). |
| preconditioner | shifted auxiliary-operator preconditioner (`HDivSchurPC`) with pressure elimination and Schur-complement fieldsplit (exact pressure-mass solve). On the velocity-buoyancy block: ONE full-multigrid cycle (V-cycles) with 6 smoothing iterations per level. The smoother is additive Schwarz on column patches, solved by LU with RCM ordering; the coarse level is solved by LU. |
| ranks | 2 MPI ranks per run; the coarse LU is PETSc's default parallel LU, MUMPS in B and C (A's installation is not recorded) |
| linear stopping test | rtol 1e-8 (for A see Section 4); at most 100 outer iterations; `ksp_converged_maxits` (a run that reaches 100 counts as "maxit") |
| nonlinear stopping test | Newton line search, at most 10 Newton steps. Eisenstat-Walker (rtol_0 = 1e-2) sets the tolerance of every linear solve (the script's `--rtol 1e-6` is overridden). KSP atol 1e-9; at most 100 iterations per linear solve; Krylov restart 30. |
| outputs | the same per-run files: `error_*.out`; for the nonlinear code also `residual_*`, `snes_error_*`, `snes_ksp_cum_*`. The same plot scripts (`plot_{dx,dz,ar,shift}.py`, called at the end of each serial command list), whose PNGs sit in each data folder. |

## 4. What is different

| item | A. Paper (May 2026) | B. Rerun (2026-10-16) | C. This test |
|---|---|---|---|
| patch | star | star and single column | star and single column |
| level smoother (6 iterations) | Chebyshev in both codes | linear: Chebyshev; nonlinear: GMRES | **GMRES in both codes** |
| outer Krylov method | linear: GMRES, left preconditioning; nonlinear: GMRES, right | linear: GMRES, left; nonlinear: FGMRES, right | **FGMRES, right, in both codes** |
| residual in the linear stopping test | preconditioned | preconditioned | **true (unpreconditioned)** |
| nonlinear advection term (`utils.Nonlinear_velocity_Irk`) | `div(outer(u, w))`, the form before commit 2cd7d82 (2026-08-17) | corrected: `div(outer(w, u))` | corrected |
| amplitude of the initial buoyancy perturbation | linear 1; nonlinear 3e-4 | linear 3e-4 (since 1e4b3d5, 2026-06-01); nonlinear 3e-4 | as B. The linear solves stop on a relative test with `ksp_atol = 0`, so this scaling does not change the linear counts. |
| flags of the linear runs | `--richardson --reordering` | `--richardson --reordering --rtol 1e-08` | `--richardson --reordering --rtol 1e-08 --mg_smoother=gmres` (new option) |
| linear rtol | not passed. The committed default on 2026-05-19 was 1e-10. But the error histories have the same length as B's star runs (explicit 1e-8) at 197 of the 206 histories present in both; the 9 others are at dt <= 5 s. So these runs effectively used 1e-8, the default committed the next day in 7771540. | 1e-8 | 1e-8 |
| rotation flag | none in either code. In both scripts `--no_rotation` was then `store_false` (passing it would switch rotation ON), so rotation was off by default. | linear: none (`LB_time_slice.py` is unchanged in this respect); nonlinear: `--no_rotation`, which is now `store_true` | as B. The physics (f = 0) is the same in A, B and C. |
| Eisenstat-Walker | set in the script's option dict | `-firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2` on the command line (the dict entries are now commented out) | as B; same values as A |
| Newton absolute tolerance | PETSc default (1e-50) | `-firedrake_1_snes_atol 1e-9` | as B |
| solver code | working tree of 2026-05-19/20, committed as 7771540 (2026-05-20 23:51) | working tree of 2026-10-16 (its `--patch` option was committed in a9a20ee, 2026-10-19) | e5b192f + the uncommitted, default-preserving `--mg_smoother` option of `LB_time_slice.py` |
| Firedrake / PETSc | the installation of May 2026, not recorded in the data. The next commit, 1e4b3d5 (2026-06-01), is "commit for all firedrake 26 changes", i.e. an update to a newer Firedrake. | `~/developer-firedrake-2026-07/venv-firedrake-07` (recorded in `run.log`) | as B |
| driver | the serial test scripts; the file times show one sweep after another | `run_claude_paper.sh` in three parts (logs `Boussinesq_implicit/claude_small_patch/run_{all_sweeps,linear,nonlinear_star}.log`): nonlinear single column with 10 concurrent runs, linear (both patches) with 4, nonlinear star with 12. One job folder per run (`jobs/job_<md5>/`: `cmd.sh`, `run.log`, `status`). | the same job folders (`run_claude_paper.sh -n`), run as ONE pool of all 854 jobs (`run_pool.sh`: 12 concurrent runs, longest first, time-out 14400 s); `mpiexec --bind-to none` as in B |
| machine load | not recorded | 4-14 concurrent 2-rank runs (two parts overlapped) | 24-30 ranks on 16 physical cores (the MPI study's np = 1 runs shared the machine), so run times are not interpreted |
| extra runs | none | none | MPI study, noise check, np = 16 diagnosis (Section 5) |
| post-processing | the per-folder plot scripts; the paper's figures | the per-folder plot scripts; `Boussinesq_implicit/claude_small_patch/compare_paper_patches.py` | the per-folder plot scripts, and `analyse.py`, which writes `figures/` (one folder per question, each with README and CSV; true counts; K and rho1 per point; comparisons with A and B; the paper's figures redrawn in `50_paper_style/`) |

Consequences for comparisons:

- **B against C, nonlinear code:** the 440 commands are byte-identical. The nonlinear runs of C are therefore a repeat
  of B; `figures/40_repeatability_nonlinear/` uses them to measure the run-to-run noise.
- **B against C, linear code:** only `--mg_smoother=gmres` differs. The smoother, the Krylov method, the preconditioning
  side and the residual norm change together, and their effects are not separated (`figures/30_smoother_linear/`).
- **A against B, linear code:** B's star commands differ from A's only by the explicit `--rtol 1e-08` (the single
  column adds `--patch vanka`). The amplitude also differs, but the linear solve does not see it. B's linear star runs
  reproduce A: the history lengths are equal at 197 of 206 points. The 9 others are at dt <= 5 s, need 28-81
  iterations and differ by 1-7.
- **A against C, nonlinear code:** the smoother, the outer Krylov method, the advection term and the Newton atol all
  changed, so the paper's nonlinear counts and those of C differ for several reasons at once.
- **Newton atol 1e-9 (B, C):** the generator documents why it is set (docstring of `bash_paper_patch.py`). It stops
  Newton before a linear solve whose right-hand side is already below the KSP atol is started. Such a solve would end at
  iteration 0, and the `ksp_converged_rate` diagnostic would then abort the run. The generator states that this does not
  change the iteration counts.

## 5. Procedure of C, step by step

Commands are run from `Boussinesq_implicit/` unless stated otherwise. Python is the Firedrake venv python.

1. **Linear solver option.** `LB_time_slice.py` gained `--mg_smoother {gmres,chebyshev,richardson}`, which preserves the
   defaults. `--richardson --mg_smoother gmres` gives:
   - GMRES(6) on the multigrid levels;
   - fieldsplit_1 `preonly` (one cycle per outer iteration);
   - outer FGMRES with right preconditioning.

   Checked with `ksp_view`, and over all 854 run logs: no PETSc option was left unused.
2. **Job folders and serial scripts.** For each code and patch:
   `./run_claude_paper.sh -n -g Bash_Generator/bash_paper_patch.py -D data -S test_paper -C patch_gmres6 -m <code> -t "dt ar dx dz" -s <suffix> -x "--patch <p> [--extra=--mg_smoother=gmres]"`
   - suffix `_deg2_gmres6_<p>` for the linear code, `_deg2_ew_gmres6_<p>` for the nonlinear code;
   - the `--extra` argument is passed for the linear code only.
3. **Running.** From `claude_patch_gmres6/`:
   - `python order_jobs.py <the 16 data folders> > joblist_main.txt` orders the jobs longest first, using B's run times;
   - `./run_pool.sh joblist_main.txt 12 14400` runs them (2026-11-16, 15:28-19:16). All 854 runs ended with exit code 0.
4. **Per-folder plots.** `./run_claude_paper.sh -p` with the same arguments as in step 2.
5. **MPI study** (from `claude_patch_gmres6/`).
   - `python make_mpi_jobs.py` writes the job folders, `data/patch_gmres6/mpi/points.csv` and one job list per np.
     There are 9 points per code (nominal, the extremes of each sweep, one single-column stress point; table in
     `figures/60_mpi/README.md`), both patches, np = 1, 2, 4, 8, 16 with np <= Nx, and three runs at np = 2 and np = 8.
     The commands are the generator's with `--ncpus <np>`.
   - `./run_mpi_chain.sh` runs np = 1 (6 at a time, alongside the main pool), then np = 2, 4, 8, 16 (12, 7, 4, 2 at a
     time), 2026-11-16 15:31-20:46.
   - The coarse LU is PETSc's own LU at np = 1 and MUMPS at np >= 2.
6. **Noise check.** `python make_noise_jobs.py; ./run_pool.sh logs/joblist_noise.txt 4`: 2 nonlinear single-column
   points with K < 1, three runs each, at np = 1 and at np = 2 with SuperLU_DIST on the coarse level. They are compared
   with the three MUMPS runs at np = 2 of the MPI study.
7. **np = 16 diagnosis.** `python make_np16_check_jobs.py; ./run_pool.sh logs/joblist_np16_check.txt 2`. The aborted
   command was rerun three ways for each patch: with converged-reason views, with SuperLU_DIST on the coarse level, and with MUMPS given
   more working space.
8. **Analysis.** `python analyse.py`, or `--only main paper smoother repeat mpi overview`, writes `figures/`.

## 6. Where the new results are

| what | where |
|---|---|
| findings (read first) | `Boussinesq_implicit/claude_patch_gmres6/README.md` |
| this procedure | `Boussinesq_implicit/claude_patch_gmres6/PROCEDURE.md` |
| figures, one folder per question, each with its own `README.md` and CSV | `Boussinesq_implicit/claude_patch_gmres6/figures/` |
| - overview of the patch ratio | `figures/00_overview/` |
| - linear sweeps (dt, ar, dx, dz) | `figures/11_linear_dt/` ... `figures/14_linear_dz/` |
| - nonlinear sweeps | `figures/21_nonlinear_dt/` ... `figures/24_nonlinear_dz/` |
| - linear GMRES(6) against Chebyshev(6) (C against B) | `figures/30_smoother_linear/` |
| - nonlinear repeatability (C against B) | `figures/40_repeatability_nonlinear/` |
| - the paper's figs. 4-12 redrawn from C, per patch | `figures/50_paper_style/{star,single_column}/{paper_conventions,corrected}/` |
| - MPI ranks, noise check, np = 16 diagnosis | `figures/60_mpi/` (all runs in `mpi_counts.csv`) |
| raw data of the main sweeps | `data/patch_gmres6/data_*` (16 folders: `*.out`, `commands.txt`, per-run `jobs/job_<md5>/{cmd.sh,run.log,status}`, the per-folder plots) |
| raw data of the MPI study | `data/patch_gmres6/mpi/<code>_<patch>/np<NN>_r<R>/`; list of all runs with commands: `data/patch_gmres6/mpi/points.csv` |
| noise check | `data/patch_gmres6/mpi/noise_check/` |
| np = 16 diagnosis | `data/patch_gmres6/mpi/np16_check/` |
| serial command lists | `data/patch_gmres6/sweeps/test_paper_*_gmres6_{star,vanka}.sh` (16 files) |
| scripts | `Boussinesq_implicit/claude_patch_gmres6/`: `run_pool.sh`, `order_jobs.py`, `make_mpi_jobs.py`, `run_mpi_chain.sh`, `make_noise_jobs.py`, `make_np16_check_jobs.py`, `analyse.py` |
| logs | `data/patch_gmres6/logs/` (job lists, `joblist_main.txt` among them, and pool logs) |
| code change | `Boussinesq_implicit/LB_time_slice.py`: option `--mg_smoother` (now committed, default gmres since 6c4b0cf) |
| earlier data used for the comparisons | A and B in `data/paper_robustness/` (unchanged) |
| follow-up 2026-11-18 (Section 7): figures and tables | `figures/70_backsub_nonlinear/` (`README.md`, `counts.csv`, `paper_style/`) |
| - raw data, back-substitution sweeps | `data/patch_gmres6/data_nonlinear_<t>_deg2_ew_gmres6_backsub_<p>/` (8 folders) |
| - raw data, SuperLU_DIST repeats at K < 1 | `data/patch_gmres6/backsub_superlu/data_nonlinear_<t>_<p>_r<1,2>/` (16 folders) |
| - pilot | `data/patch_gmres6/backsub/pilot/` |
| - serial command lists | `data/patch_gmres6/sweeps/test_paper_nonlinear_*_backsub_*.sh` (8 files) |
| - scripts, logs | `make_backsub_jobs.py`, `analyse_backsub.py`; `data/patch_gmres6/logs/joblist_backsub.txt`, `pool_backsub.log` |

## 7. Follow-up of 2026-11-18: pressure back-substitution and deterministic coarse solver

- **Same as C:** the nonlinear commands of C, byte for byte (`make_backsub_jobs.py` calls `run_claude_paper.sh -n` with
  C's generator and arguments and only a new suffix `_deg2_ew_gmres6_backsub_<p>`): parameters, GMRES(6) smoothing,
  FGMRES right, Eisenstat-Walker, MUMPS coarse LU, 2 ranks, caps and tolerances.
- **Different from C:** the script `Nonlinear_Boussinesq_Irk_SC.py` at HEAD (commit 6c4b0cf, 2026-11-17), whose
  auxiliary-operator preconditioner applies the pressure back-substitution k_p <- k_p - Pi(div k_u)/delta by default
  (class `HDivSchurPCBacksub`; `--no_backsub` gives C's preconditioner). The same commit also redefines deltax (not used
  by the solve) and printed output. A 7-run pilot confirmed: HEAD = the 2026-10-16 `--pc_fix` script (same
  counts), `--no_backsub` = C's counts, two SuperLU_DIST repeats bitwise identical, no unused PETSc options.
- **Step 1:** the four sweeps for both patches (440 runs). **Step 2:** every point with K < 1 (106 per patch), with
  `-firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_pc_factor_mat_solver_type superlu_dist` appended, run twice (424
  runs). All 864 runs in one pool (`run_pool.sh logs/joblist_backsub.txt 12 14400`, 12 jobs x 2 ranks); the pool was
  interrupted once after 544 runs and restarted (finished runs are skipped); all 864 exit 0.
- **Analysis:** `python analyse_backsub.py` (Firedrake venv python) writes `figures/70_backsub_nonlinear/`. Results:
  `README.md`, Finding 6.

