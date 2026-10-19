#!/usr/bin/env bash
# Re-run of the paper's Table 5 (fieldsplit Schur-complement MG against the monolithic MG, both on the shifted
# operator A_delta) and the diagnostics around it.  One run at a time (timings), each in runs/<tag>/ with cmd.sh,
# run.log (converged_reason of the three solves: step 1, warm-up, official run) and log.txt (-log_view).
# Settings of the table: nominal geometry 300 km x 4 km (eps = 1.3e-2), dt = 100 s, delta = 1e-4, rtol 1e-8, RCM.
#   FS   = fieldsplit Schur (--richardson: outer gmres LEFT, inner preonly V-cycle, chebyshev(6) + star patch)
#   MONO = --monolithic with chebyshev(6) smoother, outer gmres LEFT (the 2026-06-01 configuration behind the table)
# Usage: run_table5.sh [set ...]   sets: A (table reproduction, 1/2/16 ranks), Aprime (true-residual norm, 2 ranks),
#        B (inner solves converged, 2 ranks), C (larger meshes); default: all
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
FS="--richardson"
MONO="--monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
: > "$HERE/commands.txt"
add() { echo "$1|$2|$3" >> "$HERE/commands.txt"; }   # tag | nprocs | solver args
SETS="${*:-A Aprime B C}"
for set in $SETS; do
  case $set in
    A)      for np in 2 1 16; do for m in "10 10" "20 20" "40 40" "40 80"; do set -- $m
              add "A_fs_np${np}_nx$1_nz$2"   $np "--nx $1 --nz $2 $FS"
              add "A_mono_np${np}_nx$1_nz$2" $np "--nx $1 --nz $2 $MONO"; done; done ;;
    Aprime) for m in "10 10" "20 20" "40 40" "40 80"; do set -- $m
              add "Ap_fs_right_np2_nx$1_nz$2"        2 "--nx $1 --nz $2 $FS --pc_side right"
              add "Ap_monocheb_right_np2_nx$1_nz$2"  2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --pc_side right"
              add "Ap_monogmres_right_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother gmres --pc_side right"
              add "Ap_monocheb_left_rich_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left -firedrake_2_ksp_norm_type unpreconditioned"; done ;;
    B)      for m in "20 20" "40 40"; do set -- $m; for tol in 1e-6 1e-10; do
              add "B_fs_inner${tol}_np2_nx$1_nz$2"   2 "--nx $1 --nz $2 $FS --inner_rtol $tol --pc_side right"
              add "B_mono_inner${tol}_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev --inner_rtol $tol --pc_side right"; done; done ;;
    C)      for m in "80 80" "80 160"; do set -- $m; for np in 16 2; do
              add "C_fs_np${np}_nx$1_nz$2"   $np "--nx $1 --nz $2 $FS"
              add "C_mono_np${np}_nx$1_nz$2" $np "--nx $1 --nz $2 $MONO"; done; done
            add "C_fs_np1_nx80_nz80"   1 "--nx 80 --nz 80 $FS"
            add "C_mono_np1_nx80_nz80" 1 "--nx 80 --nz 80 $MONO" ;;
  esac
done
echo "[$(date)] $(wc -l < "$HERE/commands.txt") runs, sequential"
while IFS='|' read -r tag np args; do
  d="$HERE/runs/$tag"; mkdir -p "$d"
  if [ -f "$d/done" ]; then echo "[$(date +%H:%M:%S)] $tag already done"; continue; fi
  if [ "$np" -eq 1 ]; then cmd="python $PY $args $BASE"; else cmd="mpiexec --bind-to none -n $np python $PY $args $BASE"; fi
  echo "$cmd" > "$d/cmd.sh"
  t0=$(date +%s); (cd "$d" && timeout 7200 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd"
  echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  its=$(grep -o "iterations [0-9]*" "$d/run.log" | tail -1)
  echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s ${its:-no-convergence-line}"
done < "$HERE/commands.txt"
echo "[$(date)] table5 runs finished"
