#!/usr/bin/env bash
# Run a single job directory created by run_claude_paper.sh.
#   usage: claude_paper_run_one.sh <job_dir> <data_dir> <timeout_seconds>
# Runs <job_dir>/cmd.sh inside <job_dir> (so that the .pvd/.vtu and PETSc text
# files of concurrent runs never collide), copies the *.out result files up into
# <data_dir> where the plot script expects them, and records the outcome:
#   cmd.done   + done     the command that completed with output (skipped on re-runs)
#   cmd.failed + failed   the command that crashed / timed out (skipped on re-runs;
#                         delete cmd.failed to retry).  A nonlinear run whose Newton
#                         iteration diverges raises firedrake's ConvergenceError before
#                         any .out file is written, so this is how a diverged run is
#                         told apart from one that never ran; plot_claude_paper.py reads
#                         the status file and reports it as 'diverged'.
#   status                "<rc> <seconds> <reason>", reason in ok|timeout|ConvergenceError|DIVERGED_*|Traceback|Killed
# A job is skipped only if its marker was written for the *same* command line, so a
# regenerated sweep with different parameters is never served stale results.
set -u
J="$1"; D="$2"; T="${3:-7200}"
cd "$J" || exit 1
if [ -f cmd.done ] && cmp -s cmd.sh cmd.done; then
    echo "[$(date +%H:%M:%S)] $(basename "$J") already done, skipped :: $(grep -o '#.*' cmd.sh)"
    exit 0
fi
if [ -f cmd.failed ] && cmp -s cmd.sh cmd.failed; then
    echo "[$(date +%H:%M:%S)] $(basename "$J") previously failed ($(cat status 2>/dev/null)), skipped; delete cmd.failed to retry :: $(grep -o '#.*' cmd.sh)"
    exit 0
fi
rm -f done failed cmd.done cmd.failed status ./*.out
start=$(date +%s)
timeout "$T" bash cmd.sh > run.log 2>&1
rc=$?
end=$(date +%s)
# ! Prefer the solver's own verdict over the generic traceback header, which comes first in the log.
reason=$(grep -oE 'ConvergenceError|DIVERGED_[A-Z_]+' run.log | head -1)
[ -z "$reason" ] && reason=$(grep -m1 -oE 'Traceback|Killed|MemoryError' run.log)
[ "$rc" -eq 124 ] && reason="timeout"
[ "$rc" -eq 0 ] && reason="ok"
echo "$rc $((end - start)) ${reason:-error}" > status
# ! The solver scripts write their VTK output (~5 MB per run) into the cwd; nothing reads it, so it is
# ! removed here unless KEEP_VTK=1 is exported (the disk was at 98 % when this was added, 2026-10-16).
if [ "${KEEP_VTK:-0}" != "1" ]; then
    rm -rf ./Nonlinear_slice_imp_ASM ./Nonlinear_slice_imp_ASM.pvd ./lb_slice_imp_ASM ./lb_slice_imp_ASM.pvd
fi
n_out=0
for f in *.out; do
    [ -e "$f" ] || continue
    cp -f "$f" "$D/" && n_out=$((n_out + 1))
done
if [ "$rc" -eq 0 ] && [ "$n_out" -gt 0 ]; then
    cp cmd.sh cmd.done
    touch done
else
    cp cmd.sh cmd.failed
    touch failed
fi
echo "[$(date +%H:%M:%S)] $(basename "$J") rc=$rc ${reason:-error} $((end - start))s outfiles=$n_out :: $(sed -n 's/.*--dt \([^ ]*\).*--shift \([^ ]*\).*/dt=\1 shift=\2/p' cmd.sh) $(grep -o '#.*' cmd.sh)"
