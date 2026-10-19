#!/usr/bin/env bash
# Diagnostic runs for item 3 (why the nonlinear code's FGMRES count at fixed K rises with dt while the
# linear code's does not).  Nominal geometry (20 x 20, 300 km x 4 km), star patch, f = 0, no EW,
# snes_atol 1e-9, delta = K * delta_thr(dt) with K = 1 unless stated; 2 MPI ranks; each run in its own
# directory runs/<tag>/ with the solver's .out files, run.log and the command.
# Usage: run_diag.sh [NPAR]      (default 4 concurrent runs)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
NPAR="${1:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/Nonlinear_Boussinesq_Irk_SC.py"
COMMON="--nx 20 --nz 20 --length 300000.0 --height 4000.0 --maxit 100 --dt_test --no_rotation -firedrake_1_snes_atol 1e-09"

shift_of() { python - "$1" "$2" << 'PYEOF'
import sys; sys.path.insert(0, sys.argv[0] and __import__('os').environ['REPO'])
import claude_paper_shift as cps
print(repr(cps.shift_value(float(sys.argv[2]), float(sys.argv[1]))))
PYEOF
}
export REPO
: > "$HERE/commands.txt"
add() { # tag dt K extra...
    local tag="$1" dt="$2" K="$3"; shift 3
    local sh; sh=$(shift_of "$dt" "$K")
    local tmax; tmax=$(python -c "print(repr(2*float('$dt')))")
    echo "$tag|mpiexec --bind-to none -n 2 python $PY $COMMON --dt $dt --tmax $tmax --shift $sh $*" >> "$HERE/commands.txt"
}
# D1 linear forms over dt (first Newton step = like-for-like with LB_time_slice.py)
for dt in 1.0 50.0 500.0 5000.0 10000.0 100000.0; do add "D1_linforms_dt$dt" $dt 1.0 --linear_forms; done
# D2 linear forms + Chebyshev smoother (the linear code's smoother)
for dt in 50.0 5000.0 100000.0; do add "D2_linforms_cheb_dt$dt" $dt 1.0 --linear_forms --mg_smoother chebyshev; done
# D3 linear forms + exact inner solve (LU on the velocity-buoyancy block)
for dt in 50.0 5000.0 100000.0; do add "D3_linforms_fs1lu_dt$dt" $dt 1.0 --linear_forms --fs1_lu; done
# D4 linear forms + outer gmres, left preconditioning (the linear code's outer solver)
for dt in 50.0 5000.0 100000.0; do add "D4_linforms_gmresleft_dt$dt" $dt 1.0 --linear_forms --outer_gmres_left; done
# D5 linear forms, smaller K (K_eff = 2K because of the stage scaling, see the write-up)
add "D5_linforms_K0.5_dt5000.0" 5000.0 0.5 --linear_forms
add "D5_linforms_K0.1_dt5000.0" 5000.0 0.1 --linear_forms
# D6 nonlinear forms + exact inner solve
add "D6_nonlin_fs1lu_dt5000.0" 5000.0 1.0 --fs1_lu
add "D6_nonlin_fs1lu_dt50.0" 50.0 1.0 --fs1_lu
# D7 nonlinear forms, amplitude of the data
add "D7_nonlin_amp3e-5_dt5000.0" 5000.0 1.0 --amplitude 3e-5
add "D7_nonlin_amp3e-3_dt5000.0" 5000.0 1.0 --amplitude 3e-3
# D8 linear forms, tighter tolerance
add "D8_linforms_rtol1e-8_dt5000.0" 5000.0 1.0 --linear_forms --rtol 1e-8
# D9 linear forms + chebyshev + gmres-left (both linear-code settings together)
for dt in 50.0 5000.0 100000.0; do add "D9_linforms_cheb_gmresleft_dt$dt" $dt 1.0 --linear_forms --mg_smoother chebyshev --outer_gmres_left; done

run_one() {
    local line="$1"; local tag="${line%%|*}"; local cmd="${line#*|}"
    local d="$HERE/runs/$tag"; mkdir -p "$d"; cd "$d" || return 1
    if [ -f done ]; then echo "[$(date +%H:%M:%S)] $tag already done"; return 0; fi
    echo "$cmd" > cmd.sh
    local t0=$(date +%s)
    bash cmd.sh > run.log 2>&1; local rc=$?
    rm -rf Nonlinear_slice_imp_ASM Nonlinear_slice_imp_ASM.pvd
    echo "$rc $(( $(date +%s) - t0 ))" > status
    [ $rc -eq 0 ] && touch done
    echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s cum=$(cat snes_ksp_cum_*.out 2>/dev/null | tr '\n' ' ')"
}
export -f run_one; export HERE
echo "[$(date)] $(wc -l < "$HERE/commands.txt") diagnostic runs, $NPAR concurrent"
xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {} < "$HERE/commands.txt"
echo "[$(date)] diagnostics finished"
