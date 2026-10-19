#!/usr/bin/env bash
# 08_solver_degree_runs.sh
# A few runs of the REAL solver (Nonlinear_Boussinesq_Irk_SC.py, unmodified) that test the
# predictions of theory_kernel.md which the production sweeps (degree 2 only) cannot test:
#   * degree 1: the single-column patch is NOT robust (K_0 ~ 4.5 c/dx^2), the star patch is;
#   * degree 2: the single-column patch saturates as K -> 0 (c -> infinity).
# Nominal configuration: Lx = 300 km, Lz = 4 km, nz = 20, dt = 50 s, f = 0; shift = K * 1.096622711232e-08.
# Each run: 2 MPI ranks, lowest priority, own directory solver_runs/<label>/, 15 min timeout.
# usage: ./08_solver_degree_runs.sh            (runs every case that has no result yet)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/mh2520/molin/phd/yr2/Boussinesq_implicit
source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH="$REPO${PYTHONPATH:+:$PYTHONPATH}"
DT0=1.096622711232e-08

run () {  # label degree patch nx K [extra args]
  local label="$1" deg="$2" patch="$3" nx="$4" K="$5"; shift 5
  local d="$HERE/solver_runs/$label"
  if [ -f "$d/result.txt" ]; then echo "[$label] already done: $(cat $d/result.txt)"; return; fi
  rm -rf "$d"; mkdir -p "$d"; cd "$d" || exit 1
  local shift_val
  shift_val=$(python -c "print(repr($K*$DT0))")
  local cmd=(mpiexec --bind-to none -n 2 python "$REPO/Nonlinear_Boussinesq_Irk_SC.py"
     --nx "$nx" --nz 20 --length 300000.0 --height 4000.0 --dt 50.0 --tmax 100.0
     --shift "$shift_val" --maxit 100 --dx_test --no_rotation --patch "$patch" --degree "$deg"
     -firedrake_1_snes_atol 1e-09 "$@")
  printf '%q ' "${cmd[@]}" > cmd.txt; echo >> cmd.txt
  local start end rc
  start=$(date +%s)
  timeout 900 nice -n 19 "${cmd[@]}" > run.log 2>&1
  rc=$?
  end=$(date +%s)
  local cum
  cum=$(cat snes_ksp_cum_*.out 2>/dev/null | awk '{printf "%d ", $1}')
  local reason
  reason=$(grep -oE 'ConvergenceError|DIVERGED_[A-Z_]+' run.log | head -1)
  echo "degree=$deg patch=$patch nx=$nx K=$K rc=$rc wall=$((end-start))s cumulative_FGMRES_its_after_Newton_steps=[ $cum] ${reason:-}" > result.txt
  rm -rf ./*.pvd ./*.vtu ./*.pvtu Nonlinear_slice_imp_ASM 2>/dev/null
  echo "[$label] $(cat result.txt)"
}

# ---- degree 1 -------------------------------------------------------------------------
run deg1_star_nx20_K1     1 star  20 1
run deg1_vanka_nx20_K1    1 vanka 20 1
run deg1_vanka_nx20_K100  1 vanka 20 100
run deg1_star_nx20_K100   1 star  20 100
run deg1_vanka_nx20_K10   1 vanka 20 10
run deg1_star_nx20_K10    1 star  20 10
run deg1_vanka_nx40_K1    1 vanka 40 1
run deg1_star_nx40_K1     1 star  40 1
run deg1_vanka_nx10_K1    1 vanka 10 1
run deg1_star_nx10_K1     1 star  10 1
# ---- degree 2: saturation as K -> 0 ------------------------------------------------------
run deg2_vanka_nx20_K0.1  2 vanka 20 0.1
run deg2_star_nx20_K0.1   2 star  20 0.1
run deg2_vanka_nx20_K0.01 2 vanka 20 0.01
run deg2_star_nx20_K0.01  2 star  20 0.01
cd "$HERE"
for f in solver_runs/*/result.txt; do cat "$f"; done > 08_solver_degree_runs.txt
