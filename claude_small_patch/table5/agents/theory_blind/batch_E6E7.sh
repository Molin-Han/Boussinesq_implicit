#!/bin/bash
# E7: mesh variation; E6: shift variation. Variants: FS, MONO (as in Table 5), MONO with -pc_use_amat 0 on monoshiftpc
cd /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_blind
CR="-firedrake_2_ksp_converged_reason"
M="--monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
NA="-firedrake_2_monoshiftpc_pc_use_amat 0"
for mesh in "4 4" "8 4" "4 8"; do
  set -- $mesh; nx=$1; nz=$2
  B="--nx $nx --nz $nz --richardson --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing"
  ./run.sh E7_FS_nx${nx}_nz${nz} $B $CR > /dev/null
  ./run.sh E7_MONO_nx${nx}_nz${nz} $B $M $CR > /dev/null
  ./run.sh E7_MONOnoamat_nx${nx}_nz${nz} $B $M $CR $NA > /dev/null
done
for sh in 1e-3 1e-5 1e-6; do
  B="--nx 8 --nz 8 --richardson --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift $sh --maxit 100 --rtol 1.0e-8 --reordering --timing"
  ./run.sh E6_FS_shift${sh} $B $CR > /dev/null
  ./run.sh E6_MONO_shift${sh} $B $M $CR > /dev/null
  ./run.sh E6_MONOnoamat_shift${sh} $B $M $CR $NA > /dev/null
done
for d in runs/E7_* runs/E6_*; do echo "$(basename $d): $(grep 'converged due\|diverged due' $d/log.txt | grep -o 'iterations [0-9]*\|DIVERGED[A-Z_]*' | tr '\n' ' ') $(tail -1 $d/log.txt)"; done
