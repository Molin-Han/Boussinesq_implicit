#!/usr/bin/env bash
# Third batch: the corrected preconditioner output (--pc_fix) in the diagnostic copy Nonlinear_pcfix.py.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"; NPAR="${1:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH" REPO HERE
PY="$HERE/Nonlinear_pcfix.py"
COMMON="--nx 20 --nz 20 --length 300000.0 --height 4000.0 --maxit 100 --dt_test --no_rotation -firedrake_1_snes_atol 1e-09"
shift_of() { python -c "import sys; sys.path.insert(0,'$REPO'); import claude_paper_shift as c; print(repr(c.shift_value($2, $1)))"; }
: > "$HERE/commands3.txt"
add() { local tag="$1" dt="$2" K="$3"; shift 3; echo "$tag|mpiexec --bind-to none -n 2 python $PY $COMMON --dt $dt --tmax $(python -c "print(repr(2*$dt))") --shift $(shift_of $dt $K) $*" >> "$HERE/commands3.txt"; }
add "F0_copy_check_dt5000.0" 5000.0 1.0                                   # the copy without --pc_fix must reproduce 15+14
for dt in 1.0 50.0 500.0 5000.0 10000.0 100000.0; do add "F1_pcfix_dt$dt" $dt 1.0 --pc_fix; done
add "F2_pcfix_linforms_dt5000.0" 5000.0 1.0 --pc_fix --linear_forms
add "F3_pcfix_fs1lu_dt5000.0" 5000.0 1.0 --pc_fix --fs1_lu
add "F4_pcfix_K10_dt5000.0" 5000.0 10.0 --pc_fix
add "F4_pcfix_K0.1_dt5000.0" 5000.0 0.1 --pc_fix
add "F5_pcfix_vanka_dt5000.0" 5000.0 1.0 --pc_fix --patch vanka
add "F5_pcfix_vanka_dt50.0" 50.0 1.0 --pc_fix --patch vanka
run_one() {
    local line="$1"; local tag="${line%%|*}"; local cmd="${line#*|}"
    local d="$HERE/runs/$tag"; mkdir -p "$d"; cd "$d" || return 1
    if [ -f done ]; then echo "[$(date +%H:%M:%S)] $tag already done"; return 0; fi
    echo "$cmd" > cmd.sh; local t0; t0=$(date +%s)
    bash cmd.sh > run.log 2>&1; local rc=$?
    rm -rf ./Nonlinear_slice_imp_ASM ./Nonlinear_slice_imp_ASM.pvd
    echo "$rc $(( $(date +%s) - t0 ))" > status; [ $rc -eq 0 ] && touch done
    echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s cum=$(cat snes_ksp_cum_*.out 2>/dev/null | tr '\n' ' ')"
}
export -f run_one
echo "[$(date)] $(wc -l < "$HERE/commands3.txt") runs, $NPAR concurrent"
xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {} < "$HERE/commands3.txt"
echo "[$(date)] batch 3 finished"
