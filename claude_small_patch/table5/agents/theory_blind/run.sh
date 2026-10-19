#!/bin/bash
# usage: run.sh <name> <args...>   (serial only; writes into runs/<name>/)
# Runs the UNMODIFIED repo script LB_time_slice.py with cwd = runs/<name>/ so that its
# side files (slice3D.txt from --timing) land there; stdout+stderr -> runs/<name>/log.txt
export OMP_NUM_THREADS=1
PY=/home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/python
SCRIPT=/home/mh2520/molin/phd/yr2/Boussinesq_implicit/LB_time_slice.py
D=/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_blind/runs/$1
shift
mkdir -p "$D"; cd "$D"
echo "$PY $SCRIPT $@" > cmd.txt
start=$(date +%s.%N)
timeout 175 $PY $SCRIPT "$@" > log.txt 2>&1
rc=$?
end=$(date +%s.%N)
echo "rc=$rc wall=$(echo "$end - $start" | bc)" >> log.txt
tail -1 log.txt
