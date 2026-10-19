#!/usr/bin/env bash
# The paper's four robustness sweeps (delta = C_1/dt, lists of Bash/bash.py and Bash/bash_nonlinear.py),
# re-run with the single-column (vanka) patch and, for a like-for-like baseline, with the star patch,
# for both codes.  Data folders (next to the paper's data_nonlinear_*_deg2_ew and data_linear_*_deg2):
#   ../data_nonlinear_<test>_deg2_ew_vanka   ../data_nonlinear_<test>_deg2_ew_star
#   ../data_linear_<test>_deg2_vanka         ../data_linear_<test>_deg2_star
# Usage: nohup claude_small_patch/run_all_sweeps.sh [NPAR] > claude_small_patch/run_all_sweeps.log 2>&1 &
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
NPAR="${1:-10}"
cd "$HERE" || exit 1
echo "[$(date)] start, NPAR=$NPAR"
./run_claude_paper.sh -j "$NPAR" -T 3600 -g Bash/bash_paper_patch.py -D data -S test_paper -m nonlinear -t "dx dz ar dt" -s _deg2_ew_vanka -x "--patch vanka"
./run_claude_paper.sh -j "$NPAR" -T 3600 -g Bash/bash_paper_patch.py -D data -S test_paper -m linear    -t "dx dz ar dt" -s _deg2_vanka    -x "--patch vanka"
./run_claude_paper.sh -j "$NPAR" -T 3600 -g Bash/bash_paper_patch.py -D data -S test_paper -m nonlinear -t "dx dz ar dt" -s _deg2_ew_star  -x "--patch star"
./run_claude_paper.sh -j "$NPAR" -T 3600 -g Bash/bash_paper_patch.py -D data -S test_paper -m linear    -t "dx dz ar dt" -s _deg2_star     -x "--patch star"
echo "[$(date)] all sweeps finished"
