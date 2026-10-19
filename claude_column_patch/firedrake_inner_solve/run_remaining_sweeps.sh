#!/usr/bin/env bash
# The sweeps of inner_solve.py that were not run before the first session was cut off.
# Five serial processes side by side; results_<tag>.csv next to this script.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE" || exit 1
./run_sweep.sh sweep_nx     --nx 10 20 40 80             --rhs random bump &
./run_sweep.sh sweep_nz     --nz 10 20 40 80             --rhs random bump &
./run_sweep.sh sweep_Lz     --Lz 12000 4000 1000 250     --rhs random bump &
./run_sweep.sh sweep_smooth --smooth_its 6 12 20 --patch vanka --rhs random bump &
./run_sweep.sh sweep_deg1_K --degree 1 --K 0.1 1 10 100  --rhs random bump &
wait
echo "INNER SWEEPS DONE $(date)"
