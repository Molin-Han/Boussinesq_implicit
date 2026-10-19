#!/usr/bin/env bash
# Patch enumeration with 2 and 4 MPI ranks on the base mesh of the nominal configuration (nx = 20).
source /home/mh2520/developer-firedrake-2026-07/venv-firedrake-07/bin/activate
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for np in 2 4; do for patch in vanka star; do
    D="$HERE/np${np}_${patch}"; mkdir -p "$D"; cd "$D" || exit 1
    timeout 1200 mpiexec --bind-to none -n $np python "$HERE/../../02_patch_dofs.py" --mode solver --patch $patch --nx 20 --nz 3 --tag np${np}_${patch} > run.log 2>&1
    echo "np=$np patch=$patch rc=$?"
done; done
echo "PATCH ENUMERATION DONE $(date)"
