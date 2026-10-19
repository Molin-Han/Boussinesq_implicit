#!/usr/bin/env bash
# Second batch of item-3 diagnostics (see run_diag.sh): outer GMRES with left preconditioning needs a smaller
# ksp_atol, otherwise the second Newton step's KSP stops at iteration 0 and the ksp_converged_rate diagnostic
# aborts the run (LAPACKgetrf INFO = 2); --atol 1e-13 only affects where the LAST Newton step's KSP stops.
# Usage: run_diag2.sh [NPAR]
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; NPAR="${1:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH" REPO HERE
PY="$REPO/Nonlinear_Boussinesq_Irk_SC.py"
COMMON="--nx 20 --nz 20 --length 300000.0 --height 4000.0 --maxit 100 --dt_test --no_rotation -firedrake_1_snes_atol 1e-09"
shift_of() { python -c "import sys; sys.path.insert(0,'$REPO'); import claude_paper_shift as c; print(repr(c.shift_value($2, $1)))"; }
: > "$HERE/commands2.txt"
add() { local tag="$1" dt="$2" K="$3"; shift 3; echo "$tag|mpiexec --bind-to none -n 2 python $PY $COMMON --dt $dt --tmax $(python -c "print(repr(2*$dt))") --shift $(shift_of $dt $K) $*" >> "$HERE/commands2.txt"; }
for dt in 50.0 5000.0 100000.0; do add "D4b_linforms_gmresleft_atol1e-13_dt$dt" $dt 1.0 --linear_forms --outer_gmres_left --atol 1e-13; done
for dt in 50.0 5000.0 100000.0; do add "D9b_linforms_cheb_gmresleft_atol1e-13_dt$dt" $dt 1.0 --linear_forms --mg_smoother chebyshev --outer_gmres_left --atol 1e-13; done
add "D10_linforms_atol1e-13_dt5000.0" 5000.0 1.0 --linear_forms --atol 1e-13                 # control: fgmres-right with the same atol
add "D11_nonlin_gmresleft_atol1e-13_dt5000.0" 5000.0 1.0 --outer_gmres_left --atol 1e-13     # nonlinear forms, gmres left
add "D11_nonlin_gmresleft_atol1e-13_dt50.0" 50.0 1.0 --outer_gmres_left --atol 1e-13
add "D12_linforms_fs1lu_gmresleft_atol1e-13_dt5000.0" 5000.0 1.0 --linear_forms --fs1_lu --outer_gmres_left --atol 1e-13   # exact inner + gmres left
run_one() {
    local line="$1"; local tag="${line%%|*}"; local cmd="${line#*|}"
    local d="$HERE/runs/$tag"; mkdir -p "$d"; cd "$d" || return 1
    if [ -f done ]; then echo "[$(date +%H:%M:%S)] $tag already done"; return 0; fi
    echo "$cmd" > cmd.sh; local t0; t0=$(date +%s)
    bash cmd.sh > run.log 2>&1; local rc=$?
    # VTK output of this run only (written by the solver into this directory)
    rm -rf ./Nonlinear_slice_imp_ASM ./Nonlinear_slice_imp_ASM.pvd
    echo "$rc $(( $(date +%s) - t0 ))" > status; [ $rc -eq 0 ] && touch done
    echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s cum=$(cat snes_ksp_cum_*.out 2>/dev/null | tr '\n' ' ')"
}
export -f run_one
echo "[$(date)] $(wc -l < "$HERE/commands2.txt") runs, $NPAR concurrent"
xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {} < "$HERE/commands2.txt"
echo "[$(date)] batch 2 finished"
