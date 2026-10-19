#!/usr/bin/env bash
# Run the solver script once, in a fresh directory runs/<label>/, at lowest CPU priority.
#
#   usage: run_case.sh <label> <nranks> <script.py> [solver arguments ...]
#
# Example (the nominal configuration of the task, vanka patch, 2 ranks):
#   ./run_case.sh vanka_np2 2 /home/mh2520/molin/phd/yr2/Boussinesq_implicit/Nonlinear_Boussinesq_Irk_SC.py \
#        --nx 20 --nz 20 --length 300000.0 --height 4000.0 --dt 50.0 --tmax 100.0 --shift 1.096622711232e-08 \
#        --maxit 100 --dx_test --no_rotation --patch vanka -firedrake_1_snes_atol 1e-09
#
# Written into runs/<label>/ : cmd.txt, run.log, ksp_view.txt, log_view.txt, *.out, walltime.txt
# The .pvd/.vtu files of the solver are deleted at the end (they are not needed here).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/mh2520/molin/phd/yr2/Boussinesq_implicit
label="$1"; np="$2"; script="$3"; shift 3

source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"     # utils.py, for the copy of the committed script

D="$HERE/runs/$label"
rm -rf "$D"; mkdir -p "$D"; cd "$D" || exit 1
cmd=(mpiexec --bind-to none -n "$np" python "$script" "$@"
     -firedrake_1_ksp_view ":ksp_view.txt:ascii_info:append"
     -firedrake_1_ksp_converged_reason -firedrake_1_snes_converged_reason
     -log_view ":log_view.txt")
printf '%q ' "${cmd[@]}" > cmd.txt; echo >> cmd.txt
start=$(date +%s.%N)
nice -n 19 "${cmd[@]}" > run.log 2>&1
rc=$?
end=$(date +%s.%N)
echo "rc=$rc wall=$(echo "$end - $start" | bc) s  load=$(cut -d' ' -f1-3 /proc/loadavg)" > walltime.txt
rm -rf ./*.pvd ./*.vtu ./*.pvtu Nonlinear_slice_imp_ASM 2>/dev/null
echo "[$label] $(cat walltime.txt)  cum_its: $(cat snes_ksp_cum_*.out 2>/dev/null | tr '\n' ' ')"
