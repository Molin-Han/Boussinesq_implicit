#!/usr/bin/env bash
# The Table-5 study after the cause was found (2026-10-19, 07:25): PCMG sets useAmat = TRUE by default
# (petsc/src/ksp/pc/impls/mg/mg.c:1999) and Firedrake's AssembledPC/AuxiliaryOperatorPC gives its inner PC the OUTER
# operator as Amat (firedrake/preconditioners/assembled.py: pc.setOperators(A, self.P.petscmat)). For --monolithic the
# outer operator is the UNSHIFTED Jacobian, so the finest-level smoother and residual of the "monolithic shifted" MG use
# the unshifted operator; only the patches and the coarse levels are shifted. -<prefix>pc_use_amat false makes the
# inner MG use the shifted matrix on every level ("monofix"). For the fieldsplit solver Amat is the (shifted)
# MatSchurComplement, consistent with the aux operator; "fspmat" uses the assembled aux operator instead (same maths).
# Sets: F (diagnostics, counts), T (the new Table 5: 16 ranks, FGMRES right, three repetitions, rep-major),
#       H (the old monolithic configuration on the new meshes, one repetition). Same layout as run_table5.sh.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/../.." && pwd)"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PATH="$HOME/developer-firedrake-2026-07/venv-firedrake-07/bin:$PATH"
PY="$REPO/LB_time_slice.py"
BASE="--length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-4 --maxit 100 --rtol 1.0e-8 --reordering --timing -firedrake_2_ksp_converged_reason -log_view :log.txt"
FS="--richardson"
MONO="--monolithic --mono_smoother chebyshev --outer_ksp gmres --pc_side left"
FIXM="-firedrake_2_monoshiftpc_pc_use_amat false"
FIXF="-firedrake_2_fieldsplit_1_helmholtzschurpc_pc_use_amat false"
RIGHT="--outer_ksp fgmres --pc_side right"
: > "$HERE/commands_fixed.txt"
add() { echo "$1|$2|$3|$4" >> "$HERE/commands_fixed.txt"; }   # tag | nprocs | solver args | timed (1 = wait until idle)
# ---- F: diagnostics (iteration counts; may overlap other jobs)
for m in "10 10" "20 20" "40 40" "40 80"; do set -- $m
  add "F_monofix_np2_nx$1_nz$2"        2 "--nx $1 --nz $2 $MONO $FIXM" 0
  add "F_monofix_right_np2_nx$1_nz$2"  2 "--nx $1 --nz $2 --monolithic --mono_smoother chebyshev $RIGHT $FIXM" 0
  add "F_fspmat_np2_nx$1_nz$2"         2 "--nx $1 --nz $2 $FS $FIXF" 0
done
for m in "20 20" "40 40"; do set -- $m
  add "F_monofix_np1_nx$1_nz$2"        1 "--nx $1 --nz $2 $MONO $FIXM" 0
  add "F_fsdirect_np1_nx$1_nz$2"       1 "--nx $1 --nz $2 $FS --direct" 0
  add "F_monofixgmres_right_np2_nx$1_nz$2" 2 "--nx $1 --nz $2 --monolithic --mono_smoother gmres $RIGHT $FIXM" 0
done
for m in "20 20" "40 40" "40 80"; do set -- $m
  add "F_monofix_np16_nx$1_nz$2"       16 "--nx $1 --nz $2 $MONO $FIXM" 0
done
# ---- T: the new Table 5 (timed), rep-major; H: old monolithic configuration, after rep 1
MESHES="20,20 40,40 80,80 80,160 160,160"
for rep in 1 2 3; do
  for m in $MESHES; do nx=${m%,*}; nz=${m#*,}
    add "T_fs_rep${rep}_np16_nx${nx}_nz${nz}"      16 "--nx $nx --nz $nz $FS $RIGHT" 1
    add "T_fspmat_rep${rep}_np16_nx${nx}_nz${nz}"  16 "--nx $nx --nz $nz $FS $RIGHT $FIXF" 1
    add "T_monofix_rep${rep}_np16_nx${nx}_nz${nz}" 16 "--nx $nx --nz $nz --monolithic --mono_smoother chebyshev $RIGHT $FIXM" 1
  done
  if [ $rep -eq 1 ]; then
    for m in "20,20" "40,40" "80,80"; do nx=${m%,*}; nz=${m#*,}
      add "H_mono_right_np16_nx${nx}_nz${nz}" 16 "--nx $nx --nz $nz --monolithic --mono_smoother chebyshev $RIGHT" 1
    done
  fi
done
wait_idle() { while pgrep -f "LB_time_slice.py" > /dev/null; do sleep 20; done; }
echo "[$(date)] $(wc -l < "$HERE/commands_fixed.txt") runs, sequential"
while IFS='|' read -r tag np args timed; do
  d="$HERE/runs/$tag"; mkdir -p "$d"
  if [ -f "$d/done" ]; then echo "[$(date +%H:%M:%S)] $tag already done"; continue; fi
  [ "$timed" = "1" ] && wait_idle
  if [ "$np" -eq 1 ]; then cmd="python $PY $args $BASE"; else cmd="mpiexec --bind-to none -n $np python $PY $args $BASE"; fi
  echo "$cmd" > "$d/cmd.sh"
  t0=$(date +%s); (cd "$d" && timeout 7200 bash cmd.sh > run.log 2>&1 < /dev/null); rc=$?
  rm -rf "$d/lb_slice_imp_ASM" "$d/lb_slice_imp_ASM.pvd"
  echo "$rc $(( $(date +%s) - t0 ))" > "$d/status"; [ $rc -eq 0 ] && touch "$d/done"
  its=$(grep "Linear firedrake_2_ solve" "$d/run.log" | grep -o "iterations [0-9]*" | awk '{printf "%s/", $2}')
  echo "[$(date +%H:%M:%S)] $tag rc=$rc $(( $(date +%s) - t0 ))s its ${its:-none}"
done < "$HERE/commands_fixed.txt"
echo "[$(date)] fixed runs finished"
