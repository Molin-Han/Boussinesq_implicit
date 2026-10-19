#!/usr/bin/env bash
# Set Y: is the monolithic patch LU accurate? PETSc's native LU does not pivot, and --reordering puts RCM on the
# saddle-point patches of the monolithic solver (Firedrake's ASMPatchPC defaults to natural ordering "to avoid zero
# pivots in saddle-point problems"). Iteration counts only (runs concurrently with other jobs; timings not used).
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
B0="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
MONO="--monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
MUMPS_M="-firedrake_2_monoshiftpc_mg_levels_pc_star_sub_sub_pc_factor_mat_solver_type mumps"
MUMPS_F="-firedrake_2_fieldsplit_1_helmholtzschurpc_mg_levels_pc_star_sub_sub_pc_factor_mat_solver_type mumps"
: > "$HERE/commands_pivot.txt"
add() { echo "$1|$2|$3" >> "$HERE/commands_pivot.txt"; }
for m in "20 20" "40 40"; do set -- $m
  add "Y_mono_natural_np2_nx$1_nz$2"        2 "--nx $1 --nz $2 $MONO"
  add "Y_mono_rcm_mumps_np2_nx$1_nz$2"      2 "--nx $1 --nz $2 $MONO --reordering $MUMPS_M"
  add "Y_mono_natural_right_np2_nx$1_nz$2"  2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --outer_ksp fgmres --pc_side right"
  add "Y_fs_rcm_mumps_np2_nx$1_nz$2"        2 "--nx $1 --nz $2 --richardson --reordering $MUMPS_F"
done
add "Y_mono_rcm_mumps_np2_nx80_nz80"        2 "--nx 80 --nz 80 $MONO --reordering $MUMPS_M"
echo "[$(date)] $(wc -l < "$HERE/commands_pivot.txt") pivot runs"
while IFS='|' read -r tag np args; do
  d="$HERE/runs/$tag"; mkdir -p "$d"
  if [ -f "$d/done" ]; then echo "[$(date +%H:%M:%S)] $tag already done"; continue; fi
  cmd="mpiexec --bind-to none -n $np python $PY $args $B0"
  echo "$cmd" > "$d/cmd.sh"
  t0=$(date +%s); (cd "$d" && timeout 7200 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd"
  echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  its=$(grep "Linear firedrake_2_ solve" "$d/run.log" | grep -o "iterations [0-9]*" | tr '\n' ' ')
  echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s ${its:-no-convergence-line}"
done < "$HERE/commands_pivot.txt"
echo "[$(date)] pivot runs finished"
