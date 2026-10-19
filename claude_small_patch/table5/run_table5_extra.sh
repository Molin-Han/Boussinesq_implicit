#!/usr/bin/env bash
# Set D (extras, 2 ranks, (20,20) and (40,40)): the Set-B pairs again with the inner KSP's converged_reason and a
# ksp_view (inner iteration counts per outer iteration; patch sizes of both solvers). Same layout as run_table5.sh.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason -log_view :log.txt -firedrake_2_ksp_view"
: > "$HERE/commands_extra.txt"
add() { echo "$1|$2|$3" >> "$HERE/commands_extra.txt"; }
for m in "20 20" "40 40"; do set -- $m
  add "D_fs_np2_nx$1_nz$2"   2 "--nx $1 --nz $2 --richardson"
  add "D_mono_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
  # left preconditioning (the table's setting) with the true residual monitored, replacing the failed Ap_*_left_rich runs
  # (GMRES does not support the unpreconditioned norm with left preconditioning)
  add "D_fs_lefttrue_np2_nx$1_nz$2"   2 "--nx $1 --nz $2 --richardson -firedrake_2_ksp_monitor_true_residual"
  add "D_mono_lefttrue_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left -firedrake_2_ksp_monitor_true_residual"
  for tol in 1e-6 1e-10; do
    add "D_fs_inner${tol}_np2_nx$1_nz$2"   2 "--nx $1 --nz $2 --richardson --inner_rtol $tol --pc_side right -firedrake_2_fieldsplit_1_ksp_converged_reason"
    add "D_mono_inner${tol}_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --inner_rtol $tol --pc_side right -firedrake_2_monoshiftpc_ksp_converged_reason"
  done
done
echo "[$(date)] $(wc -l < "$HERE/commands_extra.txt") extra runs, sequential"
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
done < "$HERE/commands_extra.txt"
echo "[$(date)] extra runs finished"
