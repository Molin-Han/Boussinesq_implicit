#!/usr/bin/env bash
# MPI-robustness study (jobs made by make_mpi_jobs.py): np=1 runs next to the main sweep pool, then,
# once the main pool has finished, np = 2, 4, 8, 16 one after the other.  Concurrency keeps the total
# number of MPI ranks at about 24-32 on the 16-core / 32-thread machine (iteration counts do not depend
# on it; timings are not measured in this study).  Longer (nonlinear) jobs are started first.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1
run() { sort -r "logs/joblist_mpi_np$1.txt" > "logs/joblist_mpi_np$1.sorted"; ./run_pool.sh "logs/joblist_mpi_np$1.sorted" "$2" 14400; }
run 01 6
while pgrep -f "run_pool.sh joblist_main.txt" > /dev/null; do sleep 120; done
run 02 12
run 04 7
run 08 4
run 16 2
echo "[$(date +%F' '%T)] MPI chain done"
