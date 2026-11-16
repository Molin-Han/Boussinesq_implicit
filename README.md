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
| `Bash/` | Command generators that print the sweeps: `bash.py` / `bash_nonlinear.py` (the paper's original linear / nonlinear sweeps), `bash_paper_patch.py` (the same sweeps with `--patch star\|vanka`), `bash_claude_paper.py` (K sweeps), `Colin_bash.py` (unrelated template) |
| `sweeps/<category>/test_*.sh` | Serial sweep scripts, each a list of `mpiexec` commands plus a plot command, sorted by the data folder they produce: `paper_robustness/` (26), `claude_paper/` (30), `early_prototypes/` (1) |
| `docs/` | Notes: `Bug_advection_operator.md`, the mountain and lock-exchange summaries; `docs/claude_studies/` has the four Claude study records (K-shift tests, column patch, small patch / pcfix, paper update and Table 5) |
| `claude_column_patch/`, `claude_small_patch/` | Claude studies with their code and results (star vs single-column patch; small patch, pcfix diagnostics, Table 5 runs, paper figure scripts) |
| `archive/` | Superseded code, kept unchanged: `old_solvers/` (LB_Irk, LB_time_slice_separate_save, Monolithic_Boussinesq, Nonlinear_Boussinesq_slice, Nonlinear_Irk_Boussinesq_direct, SLB), `small_experiments/` (mesh-checkpoint MFEs, RCM reordering), `Colin_toy/`, `Flip_MFE_test/`, `Flip_SLB_test/`, `LB_test/`, `Plot/` |

## Where the data goes

The scripts write their output into the **current folder**. By convention, every run happens
inside a data folder `yr2/data/<category>/<name>/` next to this repository (categories:
`early_prototypes`, `nonlinear_development`, `paper_robustness`, `meanflow`,
`firedrake26_checks`, `mountain`, `claude_paper`). Seen from such a folder, the code is at
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
./run_claude_paper.sh -g Bash/bash_paper_patch.py -D data -S test_paper ...                              # paper sweeps -> data/paper_robustness/
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
- `*.png`, `*.txt`, `*.out`, `*.log`, `*.vtu`, `*.pvd`, `*.h5` are git-ignored everywhere in the repository.
