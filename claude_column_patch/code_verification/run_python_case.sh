#!/usr/bin/env bash
# Run one of the verification scripts (02_patch_dofs.py, 03_operator_serial_vs_parallel.py, ...)
# at lowest CPU priority with a given number of MPI ranks.
#
#   usage: run_python_case.sh <nranks> <script.py> [arguments ...]
#
# The scripts create their own run directory below runs/; the terminal output is stored in
# runs/logs/<script>_<arguments>.log
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
np="$1"; script="$2"; shift 2
source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$HERE/runs/logs"
name="$(basename "$script" .py)_np${np}_$(echo "$*" | tr ' /' '__' | tr -d '-')"
log="$HERE/runs/logs/$name.log"
cd "$HERE" || exit 1
start=$(date +%s)
if [ "$np" -eq 1 ]; then
    nice -n 19 python "$script" "$@" > "$log" 2>&1
else
    nice -n 19 mpiexec --bind-to none -n "$np" python "$script" "$@" > "$log" 2>&1
fi
rc=$?
end=$(date +%s)
echo "[$name] rc=$rc wall=$((end - start)) s load=$(cut -d' ' -f1-3 /proc/loadavg)"
