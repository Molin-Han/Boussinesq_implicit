#!/bin/bash
# usage: run_dense.sh <name> <args...> ; serial; cwd = runs/<name>
export OMP_NUM_THREADS=1
PY=/home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/python
H=/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_blind
D=$H/runs/$1; shift; mkdir -p $D; cd $D
echo "$PY $H/spectra_driver.py $@" > cmd.txt
start=$(date +%s.%N); timeout 175 $PY $H/spectra_driver.py "$@" > log.txt 2>&1; rc=$?; end=$(date +%s.%N)
echo "rc=$rc wall=$(echo "$end - $start" | bc)" >> log.txt; tail -3 log.txt
