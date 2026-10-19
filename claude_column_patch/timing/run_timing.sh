#!/usr/bin/env bash
# Wall-clock comparison of the star patch and the single-column (vanka) patch.
#
#   usage: ./run_timing.sh [-n NRANKS] [-r REPEATS] [-o OUTDIR] [-c "CONFIG NAMES"]
#
# Runs Nonlinear_Boussinesq_Irk_SC.py --timing (no monitors) ONE RUN AT A TIME, so use it
# on a quiet machine.  The script solves the second time step twice, in the PETSc log
# stages "Warm-up-Solver" and "Second_Run"; the numbers reported are those of "Second_Run"
# (everything compiled and set up once already), read from -log_view by parse_timing.py.
# -firedrake_1_ksp_converged_reason is added so that the iteration counts of the very same
# run are in run.log.
#
# The shift is delta = K * delta_thr(dt, L_x, L_z) from claude_paper_shift.py.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
NR=2; REP=3; OUT="$HERE/runs"; ONLY=""
while getopts "n:r:o:c:" opt; do
    case $opt in
        n) NR=$OPTARG ;; r) REP=$OPTARG ;; o) OUT=$OPTARG ;; c) ONLY=$OPTARG ;;
        *) echo "usage: $0 [-n NRANKS] [-r REPEATS] [-o OUTDIR] [-c 'config names']"; exit 1 ;;
    esac
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if ! python -c "import firedrake" > /dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin/activate"
fi

# name            nx  nz  L_z      dt        K      (total FGMRES iterations star -> vanka in the sweeps)
CONFIGS=(
  "nominal        20  20  4000.0   50.0      1.0"    # 13 -> 24
  "dt1000         20  20  4000.0   1000.0    1.0"    # 20 -> 29
  "dt5000         20  20  4000.0   5000.0    1.0"    # 29 -> 31
  "dt100000       20  20  4000.0   100000.0  1.0"    # 18 -> 20
  "thin250        20  20  250.0    50.0      1.0"    #  3 -> 32
  "thin250dt1e5   20  20  250.0    100000.0  1.0"    # 40 -> 44
  "thick12000     20  20  12000.0  50.0      1.0"    # 13 -> 17
  "nx80           80  20  4000.0   50.0      1.0"    # 11 -> 15
  "nz80           20  80  4000.0   50.0      1.0"    # 11 -> 39
  "K100           20  20  4000.0   50.0      100.0"  # 62 -> 62
)
mkdir -p "$OUT"
for rep in $(seq 1 "$REP"); do
  for cfg in "${CONFIGS[@]}"; do
    read -r name nx nz lz dt K <<< "$cfg"
    if [ -n "$ONLY" ] && ! [[ " $ONLY " == *" $name "* ]]; then continue; fi
    shift_value=$(cd "$REPO" && python -c "import claude_paper_shift as c; print(repr(c.shift_value($K, $dt, 3.0e5, $lz)))")
    tmax=$(python -c "print(repr(2 * $dt))")
    for patch in star vanka; do
      D="$OUT/${name}_np${NR}_${patch}_rep${rep}"
      if [ -f "$D/finished" ]; then echo "skip $D"; continue; fi
      mkdir -p "$D"
      cmd="mpiexec --bind-to none -n $NR python $REPO/Nonlinear_Boussinesq_Irk_SC.py --nx $nx --nz $nz --length 300000.0 --height $lz --dt $dt --tmax $tmax --shift $shift_value --maxit 100 --timing --no_rotation --patch $patch -firedrake_1_snes_atol 1e-09 -firedrake_1_ksp_converged_reason -firedrake_1_snes_converged_reason -log_view :log_view.txt"
      echo "$cmd" > "$D/cmd.sh"
      start=$(date +%s.%N)
      ( cd "$D" && bash cmd.sh > run.log 2>&1 ); rc=$?
      end=$(date +%s.%N)
      echo "$rc $(python -c "print(f'{$end - $start:.2f}')")" > "$D/status"
      [ "$rc" -eq 0 ] && touch "$D/finished"
      echo "[$(date +%H:%M:%S)] $name np=$NR $patch rep=$rep rc=$rc wall=$(cut -d' ' -f2 "$D/status")s load=$(cut -d' ' -f1 /proc/loadavg)"
    done
  done
done
