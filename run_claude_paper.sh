#!/usr/bin/env bash
# Parallel driver for the parameter robustness tests with the analysis-derived shift
# delta = K * delta_thr(dt, eps, N)   (see claude_paper_shift.py).
#
#   usage: ./run_claude_paper.sh [-j NPAR] [-m "linear nonlinear"] [-t "dt ar dx dz"] [-T SECONDS] [-n] [-p] [-s SUFFIX] [-x "GENERATOR ARGS"]
#                                [-g GENERATOR] [-D DATA_PREFIX] [-S SH_PREFIX]
#     -j  concurrent runs (each uses 2 MPI ranks)           default 10
#     -m  which codes: linear (LB_time_slice.py, KSP only) and/or nonlinear (SNES)
#     -t  which sweeps: dt (shift scaling), ar, dx, dz
#     -T  per-run timeout in seconds                         default 7200
#     -n  generate the command files and job directories only, run nothing
#     -p  plot only (assumes the data is already there)
#     -s  suffix appended to the data directory / .sh names, for side experiments
#         (e.g. -s _paperK), default none
#     -x  extra arguments handed to the generator, e.g. -x "--K 1823.8 --dts 50.0"
#     -g  command generator (relative to this directory), default Bash/bash_claude_paper.py
#         (delta = K delta_thr); Bash/bash_paper_patch.py gives the paper's own sweeps
#         (delta = C_1/dt, the lists of bash.py / bash_nonlinear.py) for a chosen patch.
#     -D  prefix of the data directories, default data_claude_paper
#         -> ../<DATA_PREFIX>_<mode>_<test><SUFFIX>/   (e.g. -D data -s _deg2_ew_vanka
#         gives ../data_nonlinear_dx_deg2_ew_vanka, next to the paper's data_nonlinear_dx_deg2_ew)
#     -S  prefix of the serial .sh files written here, default test_claude_paper
#
# The VTK output every run writes into its job directory (Nonlinear_slice_imp_ASM*,
# lb_slice_imp_ASM*; ~5 MB per run, 2.2 GB for the 2026-09-29 sweeps) is deleted by
# claude_paper_run_one.sh as soon as the run has finished, unless KEEP_VTK=1 is exported.
#
# For every (mode, test) it creates ../data_claude_paper_<mode>_<test><suffix>/, writes
# the command list there (commands.txt) and as test_claude_paper_<mode>_<test><suffix>.sh
# next to this script (the serial form, same shape as test_dt.sh), runs every command in
# its own sub-directory jobs/job_<hash of the command>/ through xargs -P, collects the
# *.out files into the data directory and finally runs plot_claude_paper.py there.
#
# Job directories are named by a hash of the command line, so re-running the driver
# with a longer sweep (e.g. -x "--dts ... 200000.0 500000.0") only executes the new
# commands and the plot covers the union; finished and failed jobs are skipped (see
# claude_paper_run_one.sh).  Different K values for the same test must go to different
# suffixes, because the solver's ar/nx/nz output file names do not carry the shift.
#
# The .sh files use the ../Boussinesq_implicit/ prefix like the existing test_*.sh, so
# they can also be run serially from inside a data directory in the usual way.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
YR2="$(cd "$HERE/.." && pwd)"
NPAR=10
MODES="linear nonlinear"
TESTS="dt ar dx dz"
TIMEOUT=7200
DRY=0
PLOT_ONLY=0
SUFFIX=""
EXTRA=""
GENERATOR="Bash/bash_claude_paper.py"
DATA_PREFIX="data_claude_paper"
SH_PREFIX="test_claude_paper"
while getopts "j:m:t:T:nps:x:g:D:S:" opt; do
    case $opt in
        j) NPAR=$OPTARG ;;
        m) MODES=$OPTARG ;;
        t) TESTS=$OPTARG ;;
        T) TIMEOUT=$OPTARG ;;
        n) DRY=1 ;;
        p) PLOT_ONLY=1 ;;
        s) SUFFIX=$OPTARG ;;
        x) EXTRA=$OPTARG ;;
        g) GENERATOR=$OPTARG ;;
        D) DATA_PREFIX=$OPTARG ;;
        S) SH_PREFIX=$OPTARG ;;
        *) echo "usage: $0 [-j NPAR] [-m modes] [-t tests] [-T timeout] [-n] [-p] [-s suffix] [-x 'generator args'] [-g generator] [-D data_prefix] [-S sh_prefix]"; exit 1 ;;
    esac
done

# ! One BLAS/OpenMP thread per rank: 2 ranks x NPAR jobs already fill the cores.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# ! Use whatever `python` is on the PATH if it has Firedrake, else the local venv.
if ! python -c "import firedrake" > /dev/null 2>&1; then
    VENV="$HOME/developer-firedrake-2026-07/venv-firedrake-07"
    if [ -f "$VENV/bin/activate" ]; then
        # shellcheck disable=SC1091
        source "$VENV/bin/activate"
    else
        echo "No python with firedrake found on PATH and $VENV does not exist." >&2
        exit 1
    fi
fi
echo "python: $(command -v python)"
echo "mpiexec: $(command -v mpiexec) ($(mpiexec --version 2>&1 | head -1))"
echo "concurrency: $NPAR runs x 2 ranks on $(getconf _NPROCESSORS_ONLN) cores"

for mode in $MODES; do
    for test in $TESTS; do
        D="$YR2/${DATA_PREFIX}_${mode}_${test}${SUFFIX}"
        mkdir -p "$D/jobs"
        echo "=================================================================="
        echo "[$(date +%H:%M:%S)] $mode / $test  ->  $D   (generator $GENERATOR)"
        # shellcheck disable=SC2086
        python "$HERE/$GENERATOR" --mode "$mode" --test "$test" $EXTRA > "$D/commands.txt" \
            || { echo "generator failed"; cat "$D/commands.txt"; continue; }
        SH="$HERE/${SH_PREFIX}_${mode}_${test}${SUFFIX}.sh"
        { echo '#!/usr/bin/env bash'; echo "# generated by run_claude_paper.sh; run from inside $D like the other test_*.sh"; cat "$D/commands.txt"; } > "$SH"
        chmod +x "$SH"

        # ! One directory per run, named by a hash of the command.  Paths are made
        # ! absolute so that the command can be executed from inside its job directory,
        # ! and OpenMPI is told not to bind: its default for -n 2 is bind-to-core, which
        # ! would pin every concurrent run to cores 0 and 1.
        : > "$D/joblist.txt"
        i=0
        while IFS= read -r line; do
            case "$line" in
                mpiexec*)
                    i=$((i + 1))
                    cmd=$(printf '%s\n' "$line" \
                        | sed -e "s#\.\./Boussinesq_implicit/#$HERE/#g" \
                              -e "s#^mpiexec -n#mpiexec --bind-to none -n#")
                    h=$(printf '%s' "$cmd" | md5sum | cut -c1-12)
                    J="$D/jobs/job_$h"
                    mkdir -p "$J"
                    printf '%s\n' "$cmd" > "$J/cmd.sh"
                    echo "$J" >> "$D/joblist.txt"
                    ;;
            esac
        done < "$D/commands.txt"
        echo "$i runs generated"
        [ "$DRY" -eq 1 ] && continue

        if [ "$PLOT_ONLY" -eq 0 ]; then
            # ! Largest dt / largest K / thinnest domain first: those are the slow ones,
            # ! and starting them first keeps the tail of the sweep short.
            tac "$D/joblist.txt" \
                | xargs -P "$NPAR" -I{} bash "$HERE/claude_paper_run_one.sh" {} "$D" "$TIMEOUT" \
                | tee -a "$D/run_all.log"
            n_done=0; n_fail=0
            while IFS= read -r J; do
                [ -f "$J/cmd.done" ] && n_done=$((n_done + 1))
                [ -f "$J/cmd.failed" ] && n_fail=$((n_fail + 1))
            done < "$D/joblist.txt"
            echo "[$(date +%H:%M:%S)] $mode / $test finished: $n_done / $i runs completed with output, $n_fail failed"
        fi

        # ! Plot (the last line of commands.txt is the plot command).
        # ! (any plot_*.py: plot_claude_paper.py for the K sweeps, plot_dx.py etc. for the paper's sweeps)
        plot_cmd=$(grep -m1 '^python .*plot_[A-Za-z_]*\.py' "$D/commands.txt" | sed -e "s#\.\./Boussinesq_implicit/#$HERE/#g")
        (cd "$D" && eval "$plot_cmd") > "$D/plot.log" 2>&1 \
            && echo "plots written to $D" || { echo "PLOTTING FAILED, see $D/plot.log"; tail -5 "$D/plot.log"; }
    done
done
echo "[$(date +%H:%M:%S)] all done"
