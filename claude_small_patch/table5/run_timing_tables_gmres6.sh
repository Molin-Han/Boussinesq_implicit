#!/usr/bin/env bash
# Re-measure the paper's two timing tables with the GMRES(6) level smoother (decision of 2026-11-16: six GMRES
# iterations on the multigrid levels in every experiment). Run on the idle Linux workstation (64 GB, 16 cores),
# from any folder:  bash claude_small_patch/table5/run_timing_tables_gmres6.sh [R] [T] [B] [C]   (default: all four sets)
# Instructions for the person or Claude session running it: WORKSTATION_TIMING_RUN.md next to this script.
#
#   R  Table 2 (tab:rcm): fieldsplit solver, (N_x, N_z) = (20, 10|20|40|80), with and without RCM reordering of the
#      patch factorisations, 2 MPI ranks, 3 repetitions.
#   T  Table 3 (tab:monolithic): fieldsplit solver against the monolithic multigrid on A_delta, both with GMRES(6) +
#      star patch + RCM on the two finer levels and LU on the coarsest, FGMRES (right preconditioning), 16 MPI ranks,
#      5 meshes, 3 repetitions.
#   B  The "exact solve" sentence below Table 3: the cycle iterated to 1e-10 (FGMRES around it) at (20,20), (40,40),
#      2 MPI ranks, both solvers.
#   C  The cost of one outer iteration with the star patch against the single column (Supplementary Section
#      sm:single-column, "about 30% cheaper"): the nonlinear solver on the nominal mesh at dt = 50 s, C_1 = 1e-5
#      (stage shift 1e-7), 2 MPI ranks, 3 repetitions; the number is PCApply / (number of applications) in the
#      stage "Second_Run" of log.txt. On Molin's Mac with Firedrake 2026.11 the single column was only about 11%
#      cheaper, because the multigrid transfers dominate there; the paper's other timings come from this machine.
#
# Common settings (those of both tables): L_x = 300 km, L_z = 4 km, dt = 100 s, delta = 1e-4, rtol 1e-8 (atol 0),
# at most 100 iterations, three multigrid levels (N_x counts the columns of the coarsest level, the finest has 4 N_x).
# With --richardson the code now uses GMRES(6) on the levels and FGMRES outside (LB_time_slice.py, 2026-11-16); the
# monolithic solver uses --mono_smoother gmres and works on A_delta on every level (pc_use_amat false is the default
# since commit e5b192f).
#
# Each run goes into runs_gmres6/<tag>/ with cmd.sh, run.log (converged reasons of the three solves: step 1, warm-up,
# official run), log.txt (-log_view; PCApply of the stage "Official-Run" is the number in the tables) and status.
# Runs that finished are skipped on a second call. Summarise with parse_timing_gmres6.py (same folder).
#
# Environment overrides (all optional):
#   FIREDRAKE_VENV  the Firedrake virtualenv (default: $HOME/developer-firedrake-2026-07/venv-firedrake-07, the one
#                   the October timings used on the workstation)
#   OUT             output folder (default: runs_gmres6 next to this script)
#   ONLY            extended regular expression; only run tags that match it (e.g. ONLY='^C_' or ONLY='nz80')
#   ALLOW_BUSY      set to 1 to run although other MPI jobs run or the load is above 1.5 (timings then unusable)
#   SKIP_PREFLIGHT  set to 1 to skip all checks below (not recommended)
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
VENV="${FIREDRAKE_VENV:-$HOME/developer-firedrake-2026-07/venv-firedrake-07}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$VENV/bin:$PATH"
PY="$REPO/LB_time_slice.py"
OUT="${OUT:-$HERE/runs_gmres6}"; mkdir -p "$OUT"
ONLY="${ONLY:-.}"
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
FS="--richardson --mg_smoother gmres"
MONO="--monolithic --mono_smoother gmres --outer_ksp fgmres --pc_side right"
NLPY="$REPO/Nonlinear_Boussinesq_Irk_SC.py"
NLBASE="--nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 1e-07 --maxit 100 --no_rotation --ew --snes_atol 1e-9 --timing -firedrake_1_ksp_converged_reason -log_view :log.txt"

# ---- preflight: the right code, a working Firedrake, an idle machine --------------------------------------------
if [ "${SKIP_PREFLIGHT:-0}" != 1 ]; then
  fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 2; }
  [ -x "$VENV/bin/python" ] || fail "no python in FIREDRAKE_VENV=$VENV (set FIREDRAKE_VENV to the Firedrake virtualenv)"
  grep -q -- "--mg_smoother" "$PY" || fail "$PY has no --mg_smoother option: the code of 2026-11-16 is not here (git pull / copy it first)"
  grep -q "HDivSchurPCBacksub" "$NLPY" || fail "$NLPY has no HDivSchurPCBacksub: the code of 2026-11-16 is not here"
  grep -q -- "--snes_atol" "$NLPY" || fail "$NLPY has no --snes_atol option: the code of 2026-11-16 is not here"
  "$VENV/bin/python" -c "import firedrake, irksome" > /dev/null 2>&1 || fail "'import firedrake, irksome' fails with $VENV/bin/python"
  for script in "$PY" "$NLPY"; do      # the scripts must at least start (e.g. no conflicting command-line options)
    "$VENV/bin/python" "$script" --help > /dev/null 2>&1 || fail "'python $script --help' fails: run it by hand to see the error, and report it to Molin (do not edit the solver code)"
  done
  command -v mpiexec > /dev/null || fail "no mpiexec on PATH"
  busy=$(pgrep -f "mpiexec|prterun|mpirun" 2>/dev/null | wc -l | tr -d ' ')
  load=$(cut -d' ' -f1 /proc/loadavg 2>/dev/null || uptime | sed -E 's/.*load averages?: *([0-9.]+).*/\1/')
  echo "preflight: venv $VENV, $(nproc 2>/dev/null || sysctl -n hw.ncpu) cores, 1-min load $load, other MPI jobs: $busy"
  if [ "${ALLOW_BUSY:-0}" != 1 ] && { [ "$busy" -gt 0 ] || awk -v l="$load" 'BEGIN{exit !(l > 1.5)}'; }; then
    fail "the machine is busy (load $load, $busy MPI jobs): timings need it idle. Wait until it is idle (or ALLOW_BUSY=1 to run anyway)."
  fi
fi

: > "$OUT/commands.txt"
add() { [[ "$1" =~ $ONLY ]] && echo "$1|$2|$3" >> "$OUT/commands.txt"; return 0; }   # tag | nprocs | solver args ("NL ..." = nonlinear solver)
SETS="${*:-R T B C}"
for set in $SETS; do
  case $set in
    R) for rep in 1 2 3; do for nz in 10 20 40 80; do
         add "R_rcm_rep${rep}_np2_nx20_nz${nz}"   2 "--nx 20 --nz $nz $FS --reordering"
         add "R_norcm_rep${rep}_np2_nx20_nz${nz}" 2 "--nx 20 --nz $nz $FS"; done; done ;;
    T) for rep in 1 2 3; do for m in "20 20" "40 40" "80 80" "80 160" "160 160"; do set -- $m
         add "T_fs_rep${rep}_np16_nx$1_nz$2"   16 "--nx $1 --nz $2 $FS --reordering"
         add "T_mono_rep${rep}_np16_nx$1_nz$2" 16 "--nx $1 --nz $2 $MONO --reordering"; done; done ;;
    B) for m in "20 20" "40 40"; do set -- $m
         add "B_fs_inner1e-10_np2_nx$1_nz$2"   2 "--nx $1 --nz $2 $FS --reordering --inner_rtol 1e-10 --pc_side right"
         add "B_mono_inner1e-10_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 $MONO --reordering --inner_rtol 1e-10"; done ;;
    C) for rep in 1 2 3; do for p in star vanka; do
         add "C_${p}_rep${rep}_np2_nx20_nz20" 2 "NL --patch $p"; done; done ;;
    *) echo "unknown set '$set' (use R, T, B, C)" >&2; exit 2 ;;
  esac
done
echo "[$(date)] $(wc -l < "$OUT/commands.txt") runs, one at a time (timings: keep the machine otherwise idle)"
while IFS='|' read -r tag np args; do
  d="$OUT/$tag"; mkdir -p "$d"
  if [ -f "$d/done" ]; then echo "[$(date +%H:%M:%S)] $tag already done"; continue; fi
  if [ "${args%% *}" = "NL" ]; then
    echo "mpiexec --bind-to none -n $np python $NLPY ${args#NL } $NLBASE" > "$d/cmd.sh"
  else
    echo "mpiexec --bind-to none -n $np python $PY $args $BASE" > "$d/cmd.sh"
  fi
  t0=$(date +%s); (cd "$d" && timeout 7200 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd" "$d/Nonlinear_slice_imp_ASM" "$d/Nonlinear_slice_imp_ASM.pvd"
  echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  its=$(grep -E "Linear firedrake_[12]_ solve" "$d/run.log" | grep -o "iterations [0-9]*" | tail -1)
  echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s ${its:-no-convergence-line}"
done < "$OUT/commands.txt"
echo "[$(date)] timing runs finished; results in $OUT"
