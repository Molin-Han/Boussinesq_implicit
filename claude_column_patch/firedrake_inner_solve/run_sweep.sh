#!/usr/bin/env bash
# Wrapper used for every run of inner_solve.py:  ./run_sweep.sh <tag> [inner_solve.py options]
# Serial, one thread, never writes byte code into the repository, log in logs/<tag>.log,
# results in results_<tag>.csv (+ results_<tag>_hist.jsonl).  Refuses to overwrite results.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TAG="$1"; shift
source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$HERE/logs" "$HERE/scratch/run_$TAG"
if [ -f "$HERE/results_$TAG.csv" ]; then echo "results_$TAG.csv exists, refusing to overwrite"; exit 1; fi
cd "$HERE/scratch/run_$TAG" || exit 1
NP="${NP:-1}"
start=$(date +%s)
if [ "$NP" -gt 1 ]; then
    timeout 1500 mpiexec --bind-to none -n "$NP" python "$HERE/inner_solve.py" --tag "$TAG" --out "$HERE/results_$TAG.csv" "$@" > "$HERE/logs/$TAG.log" 2>&1
else
    timeout 1500 python "$HERE/inner_solve.py" --tag "$TAG" --out "$HERE/results_$TAG.csv" "$@" > "$HERE/logs/$TAG.log" 2>&1
fi
rc=$?
end=$(date +%s)
echo "run_sweep $TAG rc=$rc $((end - start))s" >> "$HERE/logs/$TAG.log"
echo "run_sweep $TAG rc=$rc $((end - start))s"
