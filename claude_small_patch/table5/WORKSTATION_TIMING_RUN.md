# Timing runs on the workstation (GMRES(6) smoother)

Written 2026-11-16 on Molin's Mac, for whoever runs this on the Linux workstation: Molin, or a Claude session
started there. Read all of it once before starting. The whole job is one script, about 2–3 hours, unattended.

## Why

- **The change.** On 2026-11-16 Molin decided that every experiment in the paper uses six GMRES iterations as
  the multigrid level smoother. Before that, the linear and nonlinear sweeps used six Chebyshev iterations.
- **Already redone, on the Mac:**
  - all iteration sweeps (854 runs);
  - every figure and iteration count in the paper.
- **Still from October, with the old smoother:**
  - Table 2 (`tab:rcm`, the RCM reordering);
  - Table 3 (`tab:monolithic`, fieldsplit against monolithic);
  - the "exact solve" sentence under Table 3;
  - the single column's "about 30% cheaper per outer iteration" (supplement SM3).
- **Why here.** These are timings, and the paper's timings come from this machine (2 × Xeon Silver 4215,
  16 physical cores, 62 GiB). The Mac cannot reproduce them: it has only 12 cores, and its single-column ratio
  came out at 11% in one session and 28% in another.

## Before you start

1. **The machine must be idle.**
   - Molin had a test of his own running on 2026-11-16. Wait until it has finished, and do not stop it.
   - Check with `uptime` (1-minute load well below 1) and `pgrep -fa 'mpiexec|mpirun|prterun|python'`
     (none of his jobs).
   - The script refuses to start while the load is above 1.5 or another MPI job runs.
2. **The code of 2026-11-16 must be here.** Go to this repository on the workstation, e.g.
   `/home/mh2520/molin/phd/yr2/Boussinesq_implicit`.
   - Run `git pull` if Molin has pushed.
   - The script's preflight checks for:
     - `--mg_smoother` in `LB_time_slice.py`;
     - `HDivSchurPCBacksub` and `--snes_atol` in `Nonlinear_Boussinesq_Irk_SC.py`.
   - If they are missing, stop and ask Molin to commit and push from the Mac. Do not edit the solver code
     here. The files needed are:
     - `LB_time_slice.py`, `Nonlinear_Boussinesq_Irk_SC.py`, `utils.py`;
     - `claude_small_patch/table5/run_timing_tables_gmres6.sh`, `parse_timing_gmres6.py`, and this file.
3. **The Firedrake environment.**
   - The default is `~/developer-firedrake-2026-07/venv-firedrake-07`, the one the October timings used.
   - If it is missing, find the venv (`ls ~ | grep -i firedrake`) and pass `FIREDRAKE_VENV=/path/to/venv`.
   - Irksome must import in it. The script checks `import firedrake, irksome`.
   - Note the versions for the record: `"$VENV"/bin/python -c "import petsc4py; print(petsc4py.__version__)"`.
   - On this older Firedrake, `ASMStarPC` already builds whole-column star patches by default.
     `pc_star_column 1` only keeps that behaviour in Firedrake 2026.11, so PETSc may list it as an unused
     option. That is harmless, and the same holds for `pc_vanka_column`.
4. **A two-minute check that everything works:**
   ```bash
   ONLY='^R_rcm_rep1_np2_nx20_nz10$' bash claude_small_patch/table5/run_timing_tables_gmres6.sh R
   ```
   - It should print `rc=0 ... iterations 11`. The same run on the Mac gave 11 (step 1: 13).
   - `claude_small_patch/table5/runs_gmres6/R_rcm_rep1_np2_nx20_nz10/run.log` should end with three
     `Linear firedrake_2_ solve converged` lines.
   - Leave the folder: it counts as one of the 64 runs and is skipped later.

## Run

From the repository root:

```bash
nohup bash claude_small_patch/table5/run_timing_tables_gmres6.sh > claude_small_patch/table5/run_timing_gmres6.log 2>&1 &
```

- **The 64 runs**, one at a time:

  | Set | Runs | What |
  |---|---|---|
  | R | 24 | `tab:rcm`: 2 ranks, N_z = 10–80, RCM and no RCM, 3 repetitions each |
  | T | 30 | `tab:monolithic`: 16 ranks, five meshes up to 2.05 M unknowns, fieldsplit and monolithic, 3 repetitions each |
  | B | 4 | the cycle iterated to 1e-10 at (20,20) and (40,40) |
  | C | 6 | the nonlinear solver, star against single column, 3 repetitions |

- **Time:** about 2–3 hours in all. The longest runs are the no-RCM run at N_z = 80 and the (160,160) runs:
  each takes up to about 10 minutes.
- **Progress:** one line per run in `run_timing_gmres6.log`, e.g.
  `[14:02:11] T_fs_rep1_np16_nx80_nz80 rc=0 140s iterations 11`.
- **While it runs**, keep the machine otherwise idle: no builds, no other runs. A Claude session should
  check the log now and then, not poll it constantly.
- **Interrupted?** Start the same command again. Finished runs (folders with a `done` file) are skipped.
- **A run failed** (`rc` not 0):
  - Read its `run.log`.
  - To repeat just that run, delete its folder and call the script with `ONLY='^<tag>$'` and its set letter.
  - `rc=124` means the 2-hour timeout per run was hit; that should not happen.

## Sanity checks while or after it runs

Each run's outer iteration count must be the same in all three repetitions.

| Runs | Expected |
|---|---|
| `T_fs` at (20,20) | 11 iterations, the Mac's count at the same settings (nominal mesh, dt = 100 s, C_1 = 1e-2). A difference of more than one means the solver here differs from the Mac's: stop and report. |
| `B_fs` at (20,20) | 11 on the Mac. |
| `C_star` / `C_vanka` | 8 / 10 applications in the timed step (`Second_Run`), as in the paper's sweeps on the nominal mesh. |
| PCApply spread over 3 repetitions | In October at most 9.7% (median 3.7%). Much larger means the machine was busy. |

## After the runs

1. Summarise:
   ```bash
   python3 claude_small_patch/table5/parse_timing_gmres6.py
   ```
   It writes `claude_small_patch/table5/results_gmres6.md` and `results_gmres6.csv`. The markdown contains:
   - a **Problems** section, which should say "none";
   - medians [min, max] per configuration;
   - draft LaTeX rows for both tables;
   - the affine fit of PCApply against N_z for `tab:rcm`;
   - the exact-solve counts;
   - the single column's cost per application.
2. Do not change the paper from this file alone. The paper work is described in
   `Patch-Preconditioner-for-Boussinesq/workstation_handoff_2026-11-16.md`; open that next if you are
   continuing with the paper.
3. Do not commit or push unless Molin asks.
   - The run folders contain logs only, a few MB, plus `runs_gmres6/commands.txt`.
   - If he asks, commit `runs_gmres6/`, `results_gmres6.*` and `run_timing_gmres6.log` together.

## What the numbers mean

- **PCApply:** the wall-clock time of all preconditioner applications in the timed solve, maximum over the
  ranks, including the numerical refactorisation of the patch and coarse matrices. It is the
  `PCApply` line of the `-log_view` stage `Official-Run` (`Second_Run` for set C) in each run's `log.txt`.
  - With FGMRES and right preconditioning there is one application per outer iteration.
  - Each linear run makes three solves (step 1, a warm-up, then the timed official solve from the same initial
    guess).
  - `its` is the official solve's count.
- **Columns that do not change:** Table 2's bandwidths (48 with RCM; 531, 1051, 1571, 2091 without), and
  Table 3's unknowns and patch sizes (30 N_z + 8 fieldsplit, 38 N_z + 8 monolithic). These are properties of
  the matrices, not of the smoother.
- **Settings of both tables** (unchanged from October):
  - L_x = 300 km, L_z = 4 km, dt = 100 s, delta = 1e-4 (C_1 = 1e-2);
  - FGMRES with right preconditioning to a relative residual of 1e-8;
  - three multigrid levels: N_x counts the columns of the coarsest level, the finest has 4 N_x.
