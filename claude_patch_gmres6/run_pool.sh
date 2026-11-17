#!/usr/bin/env bash
# Run a list of job directories through ONE xargs pool (better packing than one run_claude_paper.sh per sweep).
#   usage: run_pool.sh <joblist> <NPAR> [TIMEOUT_SECONDS]
# <joblist>: one absolute job directory per line, each made by run_claude_paper.sh -n (or make_mpi_jobs.py),
# i.e. <data dir>/jobs/job_<md5>/cmd.sh.  Every job is executed by ../claude_paper_run_one.sh exactly as the
# driver does it (cwd = job dir, *.out copied to the data dir, status / cmd.done / cmd.failed recorded, VTK
# deleted, finished or failed jobs skipped on a re-run).  Plotting is done afterwards per data dir
# (run_claude_paper.sh -p ... or analyse.py).
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
LIST=$1
NPAR=$2
T=${3:-14400}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if ! python -c "import firedrake" > /dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin/activate"
fi
echo "[$(date +%F' '%T)] pool start: $(wc -l < "$LIST") jobs from $LIST, $NPAR at a time, timeout ${T}s, python $(command -v python)"
xargs -P "$NPAR" -I{} bash -c 'J="$1"; D="$(dirname "$(dirname "$J")")"; bash "$2/claude_paper_run_one.sh" "$J" "$D" "$3"' _ {} "$REPO" "$T" < "$LIST"
echo "[$(date +%F' '%T)] pool done"
