#!/usr/bin/env bash
# Sequence D: patch enumeration in parallel (small mesh) and operator-level comparison
# (small and nominal mesh).  One job after the other, never more than 2 ranks.
cd "$(dirname "${BASH_SOURCE[0]}")"
while [ ! -f runs/star_np1_repeat/walltime.txt ]; do sleep 5; done
./run_python_case.sh 2 02_patch_dofs.py --mode solver --patch vanka --nx 4 --nz 3
./run_python_case.sh 2 02_patch_dofs.py --mode solver --patch star --nx 4 --nz 3
./run_python_case.sh 2 03_operator_serial_vs_parallel.py --patch vanka --nx 4 --nz 3 --tag vanka_small_np2
./run_python_case.sh 1 03_operator_serial_vs_parallel.py --patch vanka --nx 20 --nz 20 --tag vanka_nom_np1
./run_python_case.sh 2 03_operator_serial_vs_parallel.py --patch vanka --nx 20 --nz 20 --tag vanka_nom_np2_a
./run_python_case.sh 2 03_operator_serial_vs_parallel.py --patch vanka --nx 20 --nz 20 --tag vanka_nom_np2_b
./run_python_case.sh 1 03_operator_serial_vs_parallel.py --patch star --nx 20 --nz 20 --tag star_nom_np1
./run_python_case.sh 2 03_operator_serial_vs_parallel.py --patch star --nx 20 --nz 20 --tag star_nom_np2_a
./run_python_case.sh 2 03_operator_serial_vs_parallel.py --patch star --nx 20 --nz 20 --tag star_nom_np2_b
