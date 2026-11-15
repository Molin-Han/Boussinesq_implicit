#!/usr/bin/env bash
# Set Z (2026-10-19, after the fix was made permanent in LB_time_slice.py): check the new default of the monolithic
# solver without any command-line option, the --mono_use_amat switch back to the old behaviour, the --inner_rtol
# branch, the script's default monolithic configuration, and the unchanged fieldsplit solver. (20,20), 2 ranks.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
BASE="--nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
MONO="--monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
run() { tag=$1; shift; d="$HERE/runs/$tag"; mkdir -p "$d"
  echo "mpiexec --bind-to none -n 2 python $PY $* $BASE" > "$d/cmd.sh"
  t0=$(date +%s); (cd "$d" && bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd"; echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  echo "$tag rc=$rc its $(grep 'Linear firedrake_2_ solve' "$d/run.log" | grep -o 'iterations [0-9]*' | awk '{printf "%s/", $2}')"; }
run Z_mono_np2_nx20_nz20         $MONO
run Z_mono_useamat_np2_nx20_nz20 $MONO --mono_use_amat
run Z_monodefault_np2_nx20_nz20  --monolithic
run Z_mono_inner_np2_nx20_nz20   --monolithic --mono_smoother chebyshev --inner_rtol 1e-6 --pc_side right
run Z_fs_np2_nx20_nz20           --richardson
