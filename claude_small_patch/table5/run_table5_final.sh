#!/usr/bin/env bash
# Set T: the runs behind the NEW Table 5 of the paper. Both solvers with the outer solver the paper states in its
# Section 5.1 (FGMRES, right preconditioning, so both are stopped on the same true-residual reduction rtol 1e-8),
# one multigrid cycle each with Chebyshev(6) + column-star patch + RCM, 16 MPI ranks, three repetitions per run
# (timings: median of the three Official-Run PCApply times). Same layout as run_table5.sh.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
FS="--richardson --outer_ksp fgmres --pc_side right"
MONO="--monolithic --mono_smoother chebyshev --outer_ksp fgmres --pc_side right"
NP=${NP:-16}
: > "$HERE/commands_final.txt"
add() { echo "$1|$2|$3" >> "$HERE/commands_final.txt"; }
for rep in 1 2 3; do
  for m in "20 20" "40 40" "80 80" "80 160" "160 160"; do set -- $m
    add "T_fs_rep${rep}_np${NP}_nx$1_nz$2"   $NP "--nx $1 --nz $2 $FS"
    add "T_mono_rep${rep}_np${NP}_nx$1_nz$2" $NP "--nx $1 --nz $2 $MONO"
  done
done
echo "[$(date)] $(wc -l < "$HERE/commands_final.txt") final runs, sequential"
while IFS='|' read -r tag np args; do
  d="$HERE/runs/$tag"; mkdir -p "$d"
  if [ -f "$d/done" ]; then echo "[$(date +%H:%M:%S)] $tag already done"; continue; fi
  cmd="mpiexec --bind-to none -n $np python $PY $args $BASE"
  echo "$cmd" > "$d/cmd.sh"
  t0=$(date +%s); (cd "$d" && timeout 7200 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd"
  echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  its=$(grep "Linear firedrake_2_ solve" "$d/run.log" | grep -o "iterations [0-9]*" | tail -1)
  echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s ${its:-no-convergence-line}"
done < "$HERE/commands_final.txt"
echo "[$(date)] final runs finished"
