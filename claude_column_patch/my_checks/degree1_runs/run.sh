#!/usr/bin/env bash
# Degree 1 (lowest order) with both patches on the nominal configuration, K = 1, dt = 50, 2 ranks.
source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/mh2520/molin/phd/yr2/Boussinesq_implicit
for nx in 10 20 40; do for patch in star vanka; do
  ( D="$HERE/deg1_nx${nx}_${patch}"; mkdir -p "$D"; cd "$D" || exit 1
    mpiexec --bind-to none -n 2 python $REPO/Nonlinear_Boussinesq_Irk_SC.py --degree 1 --nx $nx --nz 20 --length 300000.0 --height 4000.0 --dt 50.0 --tmax 100.0 --shift 1.096622711232e-08 --maxit 100 --dx_test --no_rotation --patch $patch -firedrake_1_snes_atol 1e-09 > run.log 2>&1
    echo "degree 1 nx=$nx $patch rc=$? cumulative: $(tr '\n' ' ' < snes_ksp_cum_dt50.0_nx${nx}.out 2>/dev/null)" ) &
done; done
wait
echo "DEGREE 1 RUNS DONE $(date)"
