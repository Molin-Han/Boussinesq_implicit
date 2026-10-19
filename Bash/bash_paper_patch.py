'''
Command generator for the robustness tests of the paper, re-run with a chosen ASM patch.

The sweeps and the shift are EXACTLY those of Bash/bash.py (linear code) and
Bash/bash_nonlinear.py (nonlinear code), which produced the paper's figures:
delta = C_1 / dt with a fixed C_1, the paper's lists of dt, heights, N_x, N_z and C_1.
Only the solver patch (--patch star | vanka) and the explicit solver flags that the
scripts at HEAD need to reproduce the paper's runs are added:

  linear     --richardson --reordering (as in bash.py), --rtol 1e-8 (the default of
             LB_time_slice.py when the paper's data were made; HEAD defaults to 1e-7),
             rotation OFF (the default of LB_time_slice.py: its --no_rotation flag is
             store_false and would switch rotation ON, so it is not passed).
  nonlinear  --no_rotation (f = 0: at the time of the paper's data the flag was
             store_false and the runs were f = 0), Eisenstat-Walker through the PETSc
             options -firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2 (as in
             data_nonlinear_*_deg2_ew, the paper's data; --no_ew switches it off), and
             -firedrake_1_snes_atol 1e-9 (stops Newton before a KSP with a right-hand
             side below ksp_atol is launched, which would abort the run through the
             ksp_converged_rate diagnostic; it does not change the iteration counts).

Prints one `mpiexec` line per run followed by the paper's plotting command
(plot_dx.py / plot_dz.py / plot_ar.py / plot_shift.py with the same arguments as the
old generators), so that

    python ../Boussinesq_implicit/Bash/bash_paper_patch.py --mode nonlinear --test dx --patch vanka

is a serial test script of the same shape as test_nonlinear_dx.sh.  run_claude_paper.sh
runs it in parallel:

    ./run_claude_paper.sh -g Bash/bash_paper_patch.py -D data -S test_paper -m nonlinear \
        -t "dx dz ar dt" -s _deg2_ew_vanka -x "--patch vanka"
'''
import os
import sys
import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description='Generate the paper\'s robustness test commands (delta = C_1/dt) for a chosen patch.',
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('--mode', choices=['linear', 'nonlinear'], required=True,
                    help='linear: LB_time_slice.py (KSP only, bash.py sweeps); nonlinear: Nonlinear_Boussinesq_Irk_SC.py (bash_nonlinear.py sweeps).')
parser.add_argument('--test', choices=['dt', 'ar', 'dx', 'dz'], required=True, help='Which parameter is swept.')
parser.add_argument('--patch', choices=['star', 'vanka'], default='star',
                    help='ASM patch of the MG smoother (passed as --patch to the solver script; nothing is added for star, '
                         'whose commands then equal the old test_*.sh lines up to the explicit flags listed in the module docstring).')
parser.add_argument('--ncpus', type=int, default=2, help='MPI ranks per run.')
parser.add_argument('--maxit', type=int, default=100, help='Max KSP iterations per linear solve.')
parser.add_argument('--dts', nargs='+', type=float, default=None, help='Override the time steps of the sweep.')
parser.add_argument('--C1', type=float, default=None, help='Override C_1 of the ar/dx/dz tests (linear default 1e-3, nonlinear 1e-5).')
parser.add_argument('--C1_list', nargs='+', type=float, default=None, help='Override the C_1 values of the dt test.')
parser.add_argument('--heights', nargs='+', type=float, default=None, help='Override the heights of the ar test.')
parser.add_argument('--nxs', nargs='+', type=int, default=None, help='Override the N_x values of the dx test.')
parser.add_argument('--nzs', nargs='+', type=int, default=None, help='Override the N_z values of the dz test.')
parser.add_argument('--rtol', type=float, default=None, help='--rtol for the solver script (linear default here: 1e-8, the paper\'s; nonlinear: not passed, script default 1e-6).')
parser.add_argument('--no_ew', action='store_true', help='(nonlinear) do not switch Eisenstat-Walker on.')
parser.add_argument('--snes_atol', type=float, default=1e-9, help='(nonlinear) -firedrake_1_snes_atol; 0 or negative = do not pass it.')
parser.add_argument('--extra', type=str, default='', help='Extra text appended verbatim to every solver command (e.g. "--linear_forms").')
parser.add_argument('--script_dir', type=str, default='../Boussinesq_implicit', help='Path prefix of the solver and plot scripts.')
args = parser.parse_args()

# ! The paper's parameter lists, copied from Bash/bash.py (linear) and Bash/bash_nonlinear.py (nonlinear).
if args.mode == 'linear':
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0]
    C1_fixed = 0.001
    C1_list = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]
    heights = [12000, 8000, 4000, 2000, 1000, 500, 250]
else:
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0]
    C1_fixed = 1e-5
    C1_list = [0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001]
    heights = [12000, 8000, 4000, 2000, 1000]
nxs = [10, 20, 40, 60, 80]
nzs = [10, 20, 40, 60, 80]
length, height, nx, nz = 3.0e5, 4000, 20, 20

if args.dts is not None:
    dts = args.dts
if args.C1 is not None:
    C1_fixed = args.C1
if args.C1_list is not None:
    C1_list = args.C1_list
if args.heights is not None:
    heights = [int(h) if float(h).is_integer() else h for h in args.heights]
if args.nxs is not None:
    nxs = args.nxs
if args.nzs is not None:
    nzs = args.nzs

if args.mode == 'linear':
    script = 'LB_time_slice.py'
    extra = '--richardson --reordering'
    rtol = 1e-8 if args.rtol is None else args.rtol
    extra += f' --rtol {rtol!r}'
else:
    script = 'Nonlinear_Boussinesq_Irk_SC.py'
    extra = '--no_rotation'
    if not args.no_ew:
        extra += ' -firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2'
    if args.snes_atol is not None and args.snes_atol > 0:
        extra += f' -firedrake_1_snes_atol {args.snes_atol!r}'
    if args.rtol is not None:
        extra += f' --rtol {args.rtol!r}'
if args.patch != 'star':
    extra += f' --patch {args.patch}'
if args.extra:
    extra += ' ' + args.extra.strip()

# ! (nx, nz, length, height, C1) members of the sweep, in the loop order of the old generators.
if args.test == 'dt':
    members = [(nx, nz, length, height, C1) for C1 in C1_list]
    loop_inner = True   # old scripts loop dt outside, C1 inside
elif args.test == 'ar':
    members = [(nx, nz, length, h, C1_fixed) for h in heights]
    loop_inner = False
elif args.test == 'dx':
    members = [(n, nz, length, height, C1_fixed) for n in nxs]
    loop_inner = False
else:
    members = [(nx, n, length, height, C1_fixed) for n in nzs]
    loop_inner = False

print(f"# {args.mode} {args.test} test of the paper (delta = C_1/dt), patch = {args.patch}"
      + (f", C_1 = {C1_fixed!r}" if args.test != 'dt' else f", C_1 in {C1_list}")
      + (", f = 0" if args.mode == 'nonlinear' else ", rotation OFF")
      + (", Eisenstat-Walker " + ("OFF" if args.no_ew else "ON") if args.mode == 'nonlinear' else '')
      + f"; solver flags: {extra}")


def line(nx_, nz_, length_, height_, dt, C1):
    options = {
        "nx": nx_,
        "nz": nz_,
        "length": length_,
        "height": height_,
        "dt": dt,
        "tmax": 2 * dt,
        "shift": np.round(C1 * dt ** (-1), decimals=16),
    }
    cmd_args = []
    for key, value in options.items():
        cmd_args += ["--" + str(key), str(value)]
    return (f"mpiexec -n {args.ncpus} python {args.script_dir}/{script} " + " ".join(cmd_args)
            + f" --maxit {args.maxit} --{args.test}_test {extra}  # C1={C1!r}")


if loop_inner:
    for dt in dts:
        for (nx_, nz_, length_, height_, C1) in members:
            print(line(nx_, nz_, length_, height_, dt, C1))
else:
    for (nx_, nz_, length_, height_, C1) in members:
        for dt in dts:
            print(line(nx_, nz_, length_, height_, dt, C1))

dts_str = " ".join(map(str, dts))
if args.test == 'dx':
    print(f"python {args.script_dir}/plot_dx.py --maxit {args.maxit} --dts {dts_str} --nxs " + " ".join(map(str, nxs))
          + f" --C1 {C1_fixed} --length {length}")
elif args.test == 'dz':
    print(f"python {args.script_dir}/plot_dz.py --maxit {args.maxit} --dts {dts_str} --nzs " + " ".join(map(str, nzs))
          + f" --C1 {C1_fixed} --height {height}")
elif args.test == 'ar':
    print(f"python {args.script_dir}/plot_ar.py --maxit {args.maxit} --dts {dts_str} --heights " + " ".join(map(str, heights))
          + f" --C1 {C1_fixed} --length {length}")
else:
    print(f"python {args.script_dir}/plot_shift.py --maxit {args.maxit} --dts {dts_str} --C1_list " + " ".join(map(str, C1_list)))
