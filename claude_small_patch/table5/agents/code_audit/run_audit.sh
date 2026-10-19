#!/usr/bin/env bash
# code_audit agent: small serial checks of the Table-5 solvers (never mpiexec; nx,nz <= 8).
# Usage: run_audit.sh TAG "solver args"   -> runs/TAG/{cmd.sh,run.log,slice3D.txt}
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO=/home/mh2520/molin/phd/yr2/Boussinesq_implicit
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/python
SCRIPT=${SCRIPT:-$REPO/LB_time_slice.py}
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason"
tag=$1; shift; d="$HERE/runs/$tag"; mkdir -p "$d"
echo "$PY $SCRIPT $* $BASE" > "$d/cmd.sh"
t0=$(date +%s); (cd "$d" && timeout 175 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
echo "$tag rc=$rc $(( $(date +%s) - t0 ))s its: $(grep -o 'iterations [0-9]*' "$d/run.log" | awk '{printf "%s/", $2}')"
