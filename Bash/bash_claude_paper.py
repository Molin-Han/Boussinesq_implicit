'''
Command generator for the parameter robustness tests with the shift chosen from the
analysis of the paper (delta = K * delta_thr, see ../claude_paper_shift.py).

Prints, exactly in the style of bash.py / bash_nonlinear.py, one `mpiexec` line per
run followed by the plotting command, so that

    python ../Boussinesq_implicit/Bash/bash_claude_paper.py --mode linear --test ar > ../Boussinesq_implicit/test_claude_paper_linear_ar.sh

reproduces a serial test script of the same shape as test_ar.sh.  The parallel
driver run_claude_paper.sh calls this generator itself.

Sweeps (defaults; every list can be overridden on the command line):

    dt   shift scaling: nominal geometry, dt x K, K in {0.1, 1, 10, 1e2, 1e3, 1e4, 1e5}
    ar   aspect ratio:  L_z in {12, 8, 4, 2, 1, 0.5, 0.25} km x dt, at K = 1
    dx   horizontal resolution: N_x in {10, 20, 40, 60, 80} x dt, at K = 1
    dz   vertical resolution:   N_z in {10, 20, 40, 60, 80} x dt, at K = 1

with dt in {1, 10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5} s.  On the nominal
geometry eps = 1/75 and N = 1e-2, so eps dt N runs from 1.3e-4 to 13 and the
cross-over dt = 1/(eps N) = 7500 s sits inside the sweep; the same is true of every
aspect ratio of the ar sweep except the thinnest two, whose cross-over
(6e4 s and 1.2e5 s) is reached or just missed at dt = 1e5 s (extend with
--dts ... 200000.0 500000.0 for those).

Rotation.  The paper's Skamarock-Klemp configuration has f = 0 (symmetric velocity
block, Theorem thm:main-sym).  LB_time_slice.py declares --no_rotation with
action='store_false', so rotation is OFF by DEFAULT there and passing the flag
would switch it ON: the linear commands therefore carry no rotation flag.
Nonlinear_Boussinesq_Irk_SC.py declares it with action='store_true' (rotation ON
by default), so the nonlinear commands pass --no_rotation explicitly.  (At the time of
the paper's nonlinear data, 2026-05-20, that script also used store_false, so the old
test_nonlinear_*.sh runs were f = 0 as well.)  --rotation switches the Coriolis term ON in
both codes, for experiments only.

Eisenstat-Walker.  The paper (appendix on Irksome) says the nonlinear solves use
the Eisenstat-Walker forcing terms, and the old data_nonlinear_*_deg2_ew runs did,
but the two lines are commented out in Nonlinear_Boussinesq_Irk_SC.py at present.
--ew switches them on from the command line through the PETSc options database
(prefix firedrake_1_, the Newton solver created by Irksome's TimeStepper in that
script) without editing the script.
'''
import os
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import claude_paper_shift as cps  # noqa: E402

parser = ArgumentParser(
    description='Generate the robustness test commands with delta = K * delta_thr(dt, eps, N).',
    formatter_class=ArgumentDefaultsHelpFormatter
)
parser.add_argument('--mode', choices=['linear', 'nonlinear'], required=True,
                    help='linear: LB_time_slice.py (KSP only); nonlinear: Nonlinear_Boussinesq_Irk_SC.py (SNES).')
parser.add_argument('--test', choices=['dt', 'ar', 'dx', 'dz'], required=True, help='Which parameter is swept.')
parser.add_argument('--ncpus', type=int, default=2, help='MPI ranks per run.')
parser.add_argument('--maxit', type=int, default=100, help='Max KSP iterations per linear solve.')
parser.add_argument('--dts', nargs='+', type=float,
                    default=[1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0, 100000.0],
                    help='Time steps swept in every test.')
parser.add_argument('--K_list', nargs='+', type=float, default=[0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0],
                    help='Safety factors K in delta = K * delta_thr, swept in the dt test.')
parser.add_argument('--K', type=float, default=1.0, help='Safety factor used in the ar, dx and dz tests.')
parser.add_argument('--heights', nargs='+', type=float, default=[12000.0, 8000.0, 4000.0, 2000.0, 1000.0, 500.0, 250.0],
                    help='Domain heights L_z (m) swept in the ar test.')
parser.add_argument('--nxs', nargs='+', type=int, default=[10, 20, 40, 60, 80], help='Coarse-mesh column counts swept in the dx test.')
parser.add_argument('--nzs', nargs='+', type=int, default=[10, 20, 40, 60, 80], help='Layer counts swept in the dz test.')
parser.add_argument('--length', type=float, default=cps.LENGTH, help='Nominal L_x (m).')
parser.add_argument('--height', type=float, default=cps.HEIGHT, help='Nominal L_z (m).')
parser.add_argument('--nx', type=int, default=20, help='Nominal coarse-mesh column count.')
parser.add_argument('--nz', type=int, default=20, help='Nominal layer count.')
parser.add_argument('--N', type=float, default=cps.N_FREQ,
                    help='Buoyancy frequency (1/s) used in the shift formula.  The solver scripts take N from '
                         'utils.buo_freq() and have no option for it, so this must stay 1e-2.')
parser.add_argument('--script_dir', type=str, default='../Boussinesq_implicit',
                    help='Path prefix of the solver and plot scripts in the printed commands.')
parser.add_argument('--rotation', action='store_true',
                    help='Switch the Coriolis term ON (f = 2 Omega sin(60 deg)) in either code.  The paper\'s '
                         'configuration, and every old run behind its figures, is f = 0 (the default here).')
parser.add_argument('--ew', action='store_true',
                    help='(nonlinear) enable Eisenstat-Walker forcing via -firedrake_1_snes_ksp_ew, as in the old '
                         'data_nonlinear_*_deg2_ew runs behind the paper\'s nonlinear figures.')
parser.add_argument('--patch', choices=['star', 'vanka'], default='star',
                    help='(nonlinear) patch of the ASM smoother, passed on as --patch to Nonlinear_Boussinesq_Irk_SC.py: '
                         'star = star of the base vertex (two columns of cells, the default and the patch of all earlier '
                         'sweeps), vanka = closure of the base cell (a single column).  Nothing is added to the command '
                         'for star, so the commands of the earlier sweeps are reproduced unchanged.')
parser.add_argument('--atol', type=float, default=None,
                    help='(nonlinear) pass --atol VALUE to the solver script (its default is 1e-9).  A Newton step whose '
                         'residual is already below ksp_atol makes the KSP stop at iteration 0, and the ksp_converged_rate '
                         'diagnostic in the script then aborts the run with "LAPACKgetrf INFO = 2" (PetscLinearRegression '
                         'on a single sample).  Thin domains (small initial residual) hit this at 1e-9; 1e-14 avoids it '
                         'without asking the KSP to reach round-off.')
parser.add_argument('--snes_atol', type=float, default=None,
                    help='(nonlinear) pass -firedrake_1_snes_atol VALUE, i.e. stop Newton once the residual is below VALUE.  '
                         'With VALUE = ksp_atol (1e-9 by default) Newton never launches a KSP whose right-hand side is '
                         'already below the KSP absolute tolerance; such a solve would stop at iteration 0 and make the '
                         'ksp_converged_rate diagnostic abort ("LAPACKgetrf INFO = 2"), and since a 0-iteration KSP '
                         'produces a zero Newton update, skipping it leaves the solution unchanged.  Preferred over --atol.')
parser.add_argument('--rtol', type=float, default=None,
                    help='pass --rtol VALUE to the solver script.  The paper\'s linear figures were made when '
                         'LB_time_slice.py defaulted to 1e-8 (commit 82cd238 of 2026-07-08 changed it to 1e-7), so '
                         'use 1e-8 for like-for-like comparisons; the nonlinear script has defaulted to 1e-6 throughout.')
parser.add_argument('--shift_override', type=float, default=None,
                    help='Use this shift (s/m^2) for every run instead of K * delta_thr.  Only for regression runs '
                         'that must reproduce an old fixed-C_1 sweep exactly (delta = C_1/dt for all heights); the '
                         'plot script keys ar/dx/dz runs by height/nx/nz, so it still finds the files, but its K, '
                         'delta and C1_equiv columns then describe the nominal K, not the shift used.')
parser.add_argument('--extra', type=str, default='',
                    help='Text appended verbatim to every solver command (e.g. "--pc_fix" for the diagnostic copy '
                         'claude_small_patch/diag_item3/Nonlinear_pcfix.py, selected with --script_dir).')
args = parser.parse_args()

if abs(args.N - 1.0e-2) > 1e-15:
    sys.exit(f"--N {args.N} differs from the N = sqrt(utils.buo_freq()) = 1e-2 hard-wired in the solver scripts; "
             "the shift would be computed for a stratification the solver does not use.")

if args.mode == 'linear' and args.patch != 'star':
    sys.exit('--patch is only implemented in Nonlinear_Boussinesq_Irk_SC.py (--mode nonlinear).')

if args.mode == 'linear':
    script = 'LB_time_slice.py'
    extra = '--richardson --reordering'          # ! as in test_dt.sh etc.; rotation is already OFF by default here.
    if args.rotation:
        extra += ' --no_rotation'                # ! store_false flag: passing it switches rotation ON in LB_time_slice.py.
else:
    script = 'Nonlinear_Boussinesq_Irk_SC.py'
    extra = '' if args.rotation else '--no_rotation'   # ! f = 0 as in the paper; ON by default in this script.
    if args.patch != 'star':
        extra += f' --patch {args.patch}'
    if args.ew:
        extra += ' -firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2'
    if args.atol is not None:
        extra += f' --atol {args.atol!r}'
    if args.snes_atol is not None:
        extra += f' -firedrake_1_snes_atol {args.snes_atol!r}'
if args.rtol is not None:
    extra += f' --rtol {args.rtol!r}'
if args.extra:
    extra += ' ' + args.extra.strip()

# ! (label, nx, nz, length, height, K) for every member of the sweep.
if args.test == 'dt':
    sweep = [(f'K{K}', args.nx, args.nz, args.length, args.height, K) for K in args.K_list]
elif args.test == 'ar':
    sweep = [(f'Lz{h}', args.nx, args.nz, args.length, h, args.K) for h in args.heights]
elif args.test == 'dx':
    sweep = [(f'nx{nx}', nx, args.nz, args.length, args.height, args.K) for nx in args.nxs]
else:
    sweep = [(f'nz{nz}', args.nx, nz, args.length, args.height, args.K) for nz in args.nzs]

print(f"# {args.mode} {args.test} test, delta = K * dt / (2 C_M^2), C_M^2 = max(C_h^2, m C_v^2), "
      f"C_h = L_x/(2 pi), C_v = L_z/pi, m = 1 + (dt N)^2/4, N = {args.N}"
      + (", rotation ON" if args.rotation else ", f = 0") + (f", patch = {args.patch}" if args.patch != 'star' else "") + (", Eisenstat-Walker ON" if args.ew else "") + (f", ksp atol = {args.atol!r}" if args.atol is not None else "") + (f", snes atol = {args.snes_atol!r}" if args.snes_atol is not None else "") + (f", rtol = {args.rtol!r}" if args.rtol is not None else "") + (f", SHIFT OVERRIDDEN to {args.shift_override!r} for every run" if args.shift_override is not None else ""))
for label, nx, nz, length, height, K in sweep:
    for dt in args.dts:
        shift = cps.shift_value(K, dt, length, height, args.N) if args.shift_override is None else float(args.shift_override)
        options = {
            'nx': nx,
            'nz': nz,
            'length': length,
            'height': height,
            'dt': dt,
            'tmax': 2 * dt,
            'shift': shift,
        }
        cmd_args = []
        for key, value in options.items():
            cmd_args += ['--' + str(key), repr(value)]
        print(f"mpiexec -n {args.ncpus} python {args.script_dir}/{script} " + ' '.join(cmd_args)
              + f" --maxit {args.maxit} --{args.test}_test {extra}"
              + f"  # {label} eps*dt*N={cps.group_eps_dt_N(dt, length, height, args.N):.3g}"
              + f" regime={cps.regime(dt, length, height, args.N)}")

plot = (f"python {args.script_dir}/plot_claude_paper.py --mode {args.mode} --test {args.test}"
        f" --maxit {args.maxit} --dts " + ' '.join(repr(dt) for dt in args.dts)
        + f" --length {args.length!r} --height {args.height!r} --nx {args.nx} --nz {args.nz} --N {args.N!r}")
if args.test == 'dt':
    plot += ' --K_list ' + ' '.join(repr(K) for K in args.K_list)
elif args.test == 'ar':
    plot += f' --K {args.K!r} --heights ' + ' '.join(repr(h) for h in args.heights)
elif args.test == 'dx':
    plot += f' --K {args.K!r} --nxs ' + ' '.join(str(n) for n in args.nxs)
else:
    plot += f' --K {args.K!r} --nzs ' + ' '.join(str(n) for n in args.nzs)
print(plot)
