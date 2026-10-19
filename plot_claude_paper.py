'''
Plotting for the parameter robustness tests with the analysis-derived shift
delta = K * delta_thr(dt, eps, N) (claude_paper_shift.py).  Run from inside a
data_claude_paper_<mode>_<test>/ directory, normally by run_claude_paper.sh; the
command line is the last line of commands.txt there.

Input files are the ones the solver scripts write with --dt_test / --ar_test /
--dx_test / --dz_test:

    error_dt{dt}_{key}{value}.out          relative solution error per KSP iteration
    residual_dt{dt}_{key}{value}.out       (nonlinear only) KSP residual per iteration
    snes_error_dt{dt}_{key}{value}.out     (nonlinear only) relative error per Newton step
    snes_ksp_cum_dt{dt}_{key}{value}.out   (nonlinear only) cumulative KSP count per Newton step

with key = shift (dt test), ar (ar test), nx (dx test), nz (dz test).

Conventions (they differ slightly from plot_shift.py / plot_ar.py, whose figures are
in the paper):
  * linear: the KSP monitor is called at iterations 0..n, so a file with n+1 entries
    is n iterations, and n is what is reported (the old scripts reported n+1).  The
    history plots show all n+1 points; the last one is the error the KSP stopped at
    (the reference is the warm-up solve, and the monitored solve's final iterate
    differs from it at round-off level), which the old scripts dropped.
  * nonlinear: the count is snes_ksp_cum[-1], the total number of KSP iterations over
    all Newton steps (the old error_shift_nonlinear.png used len(error), which is the
    total plus the number of Newton steps).  The SNES-error history drops its last
    point (the converged state, exactly the reference) as plot_shift.py did.
  * a run whose KSP hit --maxit (in any Newton step) is reported as NaN, '>maxit';
    a run that crashed (Firedrake raises ConvergenceError when Newton diverges, before
    any .out file is written) is found through jobs/*/status written by
    claude_paper_run_one.sh and reported as 'diverged' / 'timeout' / 'crash'; a run
    with no file and no job record is 'missing'.

Outputs (names follow the existing plot scripts where a counterpart exists, so that
the figures can be put next to the ones in the paper):

    design_shift.png                 delta actually used vs dt, with the two asymptotes and the cross-over
    predicted_rho01.png              predicted damping factor of the (0,1) mode vs eps dt N
    summary_iterations.csv / .md     every run: parameters, group values, iteration counts, status
    error_<T>.png                    KSP iterations (linear) vs dt, one curve per swept value
    error_<T>_vs_group.png           the same against eps dt N
    error_<T>_vs_4eps2m.png          the same against 4 eps^2 m, the exact collapse variable
    error_<T>_Robust_dt{dt}.png      (linear) error history at fixed dt, one curve per swept value
    error_Robust_K{K}.png            (linear dt test) error history at fixed K, one curve per dt
    iteration_map_<T>.png            (dt test) iteration count over the (dt, K) plane with K = 1 and
                                     K_operative(dt) = 1/min(1, 4 eps^2 m) overlaid
    snes_ksp_<T>.png                 (nonlinear) total KSP iterations over all Newton steps vs dt
    snes_ksp_<T>_vs_group.png        (nonlinear) the same against eps dt N
    snes_ksp_<T>_vs_4eps2m.png       (nonlinear) the same against 4 eps^2 m
    newton_steps_<T>.png             (nonlinear) Newton steps vs dt
    snes_error_<T>_dt{dt}.png        (nonlinear) relative error vs cumulative KSP iterations at fixed dt

with <T> = shift, AR, dx, dz; for the nonlinear dt test snes_ksp_shift.png is also
saved as error_shift_nonlinear.png, the name used in the paper.
'''
import os
import sys
import glob
import re
import csv
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import claude_paper_shift as cps  # noqa: E402

parser = ArgumentParser(
    description='Plot the robustness tests run with delta = K * delta_thr(dt, eps, N).',
    formatter_class=ArgumentDefaultsHelpFormatter
)
parser.add_argument('--mode', choices=['linear', 'nonlinear'], required=True)
parser.add_argument('--test', choices=['dt', 'ar', 'dx', 'dz'], required=True)
parser.add_argument('--maxit', type=int, default=100, help='KSP iteration cap used in the runs.')
parser.add_argument('--dts', nargs='+', type=float, required=True, help='Time steps of the sweep.')
parser.add_argument('--K_list', nargs='+', type=float, default=[1.0], help='(dt test) safety factors swept.')
parser.add_argument('--K', type=float, default=1.0, help='(ar, dx, dz tests) safety factor used.')
parser.add_argument('--heights', nargs='+', type=float, default=[], help='(ar test) L_z values swept.')
parser.add_argument('--nxs', nargs='+', type=int, default=[], help='(dx test) column counts swept.')
parser.add_argument('--nzs', nargs='+', type=int, default=[], help='(dz test) layer counts swept.')
parser.add_argument('--length', type=float, default=cps.LENGTH)
parser.add_argument('--height', type=float, default=cps.HEIGHT)
parser.add_argument('--nx', type=int, default=20)
parser.add_argument('--nz', type=int, default=20)
parser.add_argument('--N', type=float, default=cps.N_FREQ)
args = parser.parse_known_args()[0]

dts = sorted(args.dts)
mode = args.mode
test = args.test
N = args.N
FS_LEG, FS_LAB, FS_TIT = 14, 18, 16


def _sci(x):
    if x == 0:
        return r'$0$'
    exp = int(np.floor(np.log10(abs(x))))
    mantissa = x / 10.0 ** exp
    if abs(mantissa - 1.0) < 1e-9:
        return rf'$10^{{{exp}}}$'
    return rf'${mantissa:.2g}\times 10^{{{exp}}}$'


def _log_close(v, target, tol=0.2):
    if v > 0 and target > 0:
        return abs(np.log(v) - np.log(target)) <= tol
    return abs(v - target) <= tol


def _find_data_path(prefix, dt, key, target, ext='.out'):
    '''
    Pick {prefix}_dt{dt}_{key}<value>{ext} whose <value> is closest in log space to
    target, as in plot_shift.py; None if nothing is within 20 % (log).  The dt in the
    file name is the repr of the float the solver script received, which is how the
    generator printed it.
    '''
    candidates = glob.glob(f'{prefix}_dt{float(dt)!r}_{key}*{ext}')
    if not candidates:
        return None
    rx = re.compile(rf'_{re.escape(key)}([^/]+){re.escape(ext)}$')
    best, best_d = None, float('inf')
    for p in candidates:
        m = rx.search(p)
        if not m:
            continue
        try:
            v = float(m.group(1))
        except ValueError:
            continue
        d = abs(np.log(v) - np.log(target)) if (v > 0 and target > 0) else abs(v - target)
        if d < best_d:
            best, best_d = p, d
    return best if best is not None and best_d <= 0.2 else None


def _load(prefix, dt, key, target):
    path = _find_data_path(prefix, dt, key, target)
    if path is None:
        return None
    try:
        return np.atleast_1d(np.loadtxt(path))
    except (OSError, ValueError):
        return None


def _scan_failed_jobs():
    '''
    Crashed / timed-out runs, from jobs/*/status ("<rc> <seconds> <reason>") written by
    claude_paper_run_one.sh, with the parameters parsed from the job's cmd.sh.
    '''
    failed = []
    for st in glob.glob('jobs/*/status'):
        try:
            rc, _secs, reason = open(st).read().split()[:3]
        except (OSError, ValueError):
            continue
        if rc == '0':
            continue
        try:
            cmd = open(os.path.join(os.path.dirname(st), 'cmd.sh')).read()
        except OSError:
            continue

        def opt(name, cast=float):
            m = re.search(rf'--{name} (\S+)', cmd)
            return cast(m.group(1)) if m else None

        failed.append(dict(dt=opt('dt'), shift=opt('shift'), height=opt('height'), length=opt('length'),
                           nx=opt('nx', int), nz=opt('nz', int), reason=reason))
    return failed


def _failure_status(rec, mem, failed):
    for f in failed:
        if f['dt'] is None or abs(f['dt'] - rec['dt']) > 1e-9 * max(1.0, rec['dt']):
            continue
        if mem['key'] == 'shift' and not _log_close(f['shift'] or 0.0, rec['delta']):
            continue
        if mem['key'] == 'ar' and not (f['height'] == mem['height'] and f['length'] == mem['length']):
            continue
        if mem['key'] == 'nx' and f['nx'] != mem['nx']:
            continue
        if mem['key'] == 'nz' and f['nz'] != mem['nz']:
            continue
        r = f['reason']
        if r.startswith('ConvergenceError') or r.startswith('DIVERGED'):
            return 'diverged'
        if r == 'timeout':
            return 'timeout'
        return 'crash'
    return None


# ! The members of the sweep: one curve each.
members = []
if test == 'dt':
    for K in args.K_list:
        members.append(dict(legend=f'K={_sci(K)}', tag=f'K{K!r}', nx=args.nx, nz=args.nz,
                            length=args.length, height=args.height, K=K, key='shift'))
elif test == 'ar':
    for h in args.heights:
        members.append(dict(legend=f'AR={_sci(h / args.length)}', tag=f'Lz{h!r}', nx=args.nx, nz=args.nz,
                            length=args.length, height=h, K=args.K, key='ar'))
elif test == 'dx':
    for nx in args.nxs:
        members.append(dict(legend=f'$N_x$={nx}', tag=f'nx{nx}', nx=nx, nz=args.nz,
                            length=args.length, height=args.height, K=args.K, key='nx'))
else:
    for nz in args.nzs:
        members.append(dict(legend=f'$N_z$={nz}', tag=f'nz{nz}', nx=args.nx, nz=nz,
                            length=args.length, height=args.height, K=args.K, key='nz'))
if not members:
    sys.exit('nothing to plot: the sweep list for this test is empty')

T_NAME = {'dt': 'shift', 'ar': 'AR', 'dx': 'dx', 'dz': 'dz'}[test]
SWEEP_TITLE = {'dt': 'varying K', 'ar': 'varying aspect ratio', 'dx': r'varying $\Delta x$', 'dz': r'varying $\Delta z$'}[test]
failed_jobs = _scan_failed_jobs()

# ! Read everything once.  records[member index][dt] -> dict.
records = []
for mem in members:
    per_dt = {}
    for dt in dts:
        L, H, K = mem['length'], mem['height'], mem['K']
        delta = cps.shift_value(K, dt, L, H, N)
        target = {'shift': delta, 'ar': H / L, 'nx': mem['nx'], 'nz': mem['nz']}[mem['key']]
        rec = dict(dt=dt, eps=H / L, delta=delta, K=K, C1_equiv=delta * dt,
                   group=cps.group_eps_dt_N(dt, L, H, N), g4=cps.group_4eps2m(dt, L, H, N),
                   regime=cps.regime(dt, L, H, N), rho01=cps.rho_01(K, dt, L, H, N),
                   K_op=cps.K_operative(dt, L, H, N),
                   its=np.nan, total_ksp=np.nan, newton=np.nan, status='missing',
                   error=None, residual=None, snes_error=None, snes_cum=None)
        error = _load('error', dt, mem['key'], target)
        if error is not None:
            rec['error'] = error
        if mode == 'linear':
            if error is not None:
                n = len(error) - 1
                rec['its'] = np.nan if n >= args.maxit else n
                rec['status'] = 'maxit' if n >= args.maxit else 'ok'
        else:
            rec['residual'] = _load('residual', dt, mem['key'], target)
            rec['snes_error'] = _load('snes_error', dt, mem['key'], target)
            rec['snes_cum'] = _load('snes_ksp_cum', dt, mem['key'], target)
            cum = rec['snes_cum']
            if cum is not None and len(cum) >= 1:
                rec['newton'] = len(cum) - 1
                per_step = np.diff(cum) if len(cum) > 1 else np.array([cum[-1]])
                rec['its'] = float(per_step[-1])        # ! KSP iterations of the last Newton step
                if np.any(per_step >= args.maxit):
                    rec['status'] = 'maxit'
                else:
                    rec['status'] = 'ok'
                    rec['total_ksp'] = float(cum[-1])
            elif error is not None:
                rec['status'] = 'no_snes_file'
        if rec['status'] == 'missing':
            fs = _failure_status(rec, mem, failed_jobs)
            if fs is not None:
                rec['status'] = fs
        per_dt[dt] = rec
    records.append(per_dt)

# ---------------------------------------------------------------------------
# Summary table (csv + markdown), the first thing to look at.
# ---------------------------------------------------------------------------
count_key = 'its' if mode == 'linear' else 'total_ksp'
with open('summary_iterations.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mode', 'test', 'member', 'K', 'nx', 'nz', 'length', 'height', 'eps', 'dt', 'eps_dt_N',
                '4eps2m', 'regime', 'delta', 'C1_equiv', 'rho01_pred', 'K_operative',
                'ksp_its' if mode == 'linear' else 'ksp_its_last_newton_step', 'total_ksp', 'newton_steps', 'status'])
    for mem, per_dt in zip(members, records):
        for dt in dts:
            r = per_dt[dt]
            w.writerow([mode, test, mem['tag'], mem['K'], mem['nx'], mem['nz'], mem['length'], mem['height'],
                        f"{r['eps']:.6g}", dt, f"{r['group']:.6g}", f"{r['g4']:.6g}", r['regime'],
                        f"{r['delta']:.6e}", f"{r['C1_equiv']:.6e}", f"{r['rho01']:.4g}", f"{r['K_op']:.4g}",
                        r['its'], r['total_ksp'], r['newton'], r['status']])

CELL = {'missing': 'missing', 'maxit': f'>{args.maxit}', 'no_snes_file': 'no_snes',
        'diverged': 'div', 'timeout': 'timeout', 'crash': 'crash'}


def _cell(r):
    if r['status'] != 'ok':
        return CELL.get(r['status'], r['status'])
    v = r[count_key]
    return f'{int(v)}' if np.isfinite(v) else 'nan'


with open('summary_table.md', 'w') as f:
    f.write(f"# {mode} {test} test: {'KSP iterations' if mode == 'linear' else 'total KSP iterations over all Newton steps'}\n\n")
    f.write(f"delta = K * dt / (2 C_M^2), N = {N}, nominal L_x = {args.length}, L_z = {args.height}; "
            f"cross-over dt = 1/(eps N) = {cps.crossover_dt(args.length, args.height, N):.4g} s "
            f"for the nominal aspect ratio.  'div' = Newton diverged (ConvergenceError), "
            f"'>{args.maxit}' = a KSP solve hit the iteration cap.\n\n")
    f.write('| dt | eps dt N (nominal) | K_operative (nominal) | '
            + ' | '.join(m['legend'].replace('$', '') for m in members) + ' |\n')
    f.write('|---|---|---|' + '---|' * len(members) + '\n')
    for dt in dts:
        f.write(f"| {dt:g} | {cps.group_eps_dt_N(dt, args.length, args.height, N):.3g} | "
                f"{cps.K_operative(dt, args.length, args.height, N):.3g} | "
                + ' | '.join(_cell(per_dt[dt]) for per_dt in records) + ' |\n')
    if mode == 'nonlinear':
        f.write('\nNewton steps:\n\n| dt | ' + ' | '.join(m['legend'].replace('$', '') for m in members) + ' |\n')
        f.write('|---|' + '---|' * len(members) + '\n')
        for dt in dts:
            f.write(f"| {dt:g} | " + ' | '.join(
                CELL.get(per_dt[dt]['status'], '-') if not np.isfinite(per_dt[dt]['newton']) else f"{int(per_dt[dt]['newton'])}"
                for per_dt in records) + ' |\n')

# ---------------------------------------------------------------------------
# Design plots: the shift actually used, and the predicted mode damping.
# ---------------------------------------------------------------------------
dt_fine = np.logspace(np.log10(min(dts)) - 0.3, np.log10(max(dts)) + 0.3, 400)
fig, ax = plt.subplots(figsize=(10, 6))
if test == 'ar':
    for mem in members:
        ax.loglog(dt_fine, [cps.shift_threshold(d, mem['length'], mem['height'], N) for d in dt_fine],
                  label=rf"$\delta_{{thr}}$, {mem['legend']}")
        ax.axvline(cps.crossover_dt(mem['length'], mem['height'], N), ls=':', color='grey', lw=0.8)
    if args.K != 1.0:
        ax.set_title(rf'$\delta_{{thr}}$ per aspect ratio; runs use $\delta = K\,\delta_{{thr}}$ with K={_sci(args.K)}',
                     fontsize=FS_TIT)
else:
    ax.loglog(dt_fine, [cps.shift_threshold(d, args.length, args.height, N) for d in dt_fine],
              'k-', lw=2, label=r'$\delta_{thr} = \Delta t / (2 C_M^2)$')
    ax.loglog(dt_fine, cps.shift_weak_branch(dt_fine, args.length), 'b--', label=r'$2\pi^2 \Delta t / L_x^2$ (weak)')
    ax.loglog(dt_fine, cps.shift_strong_branch(dt_fine, args.height, N), 'r--',
              label=r'$2\pi^2 / (\Delta t N^2 L_z^2)$ (strong)')
    ax.loglog(dt_fine, [cps.shift_threshold_baroclinic(d, args.height, N) for d in dt_fine], 'g-.',
              label=r'$\sigma_{01}/2$ (operative, eq. C1-threshold-baroclinic)')
    ax.axvline(cps.crossover_dt(args.length, args.height, N), ls=':', color='grey',
               label=r'$\Delta t = 1/(\varepsilon N)$')
    for mem in members:
        ax.loglog(dts, [cps.shift_value(mem['K'], d, mem['length'], mem['height'], N) for d in dts],
                  'o', ms=5, label=f"used, {mem['legend']}")
ax.set_xlabel(r'$\Delta t$ (s)', fontsize=FS_LAB)
ax.set_ylabel(r'$\delta$ (s m$^{-2}$)', fontsize=FS_LAB)
if not ax.get_title():
    ax.set_title(rf'Shift parameter used in the {mode} {test} test', fontsize=FS_TIT)
ax.legend(fontsize=10, ncol=2)
ax.grid(True, which='both', alpha=0.3)
fig.savefig('design_shift.png', bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 6))
for mem in members:
    ax.semilogx([cps.group_eps_dt_N(d, mem['length'], mem['height'], N) for d in dt_fine],
                [cps.rho_01(mem['K'], d, mem['length'], mem['height'], N) for d in dt_fine], label=mem['legend'])
ax.axvline(1.0, ls=':', color='grey')
ax.set_xlabel(r'$\varepsilon\,\Delta t\,N$', fontsize=FS_LAB)
ax.set_ylabel(r'predicted $\rho_{01} = \delta/(\delta+\sigma_{01})$', fontsize=FS_LAB)
ax.set_title('Predicted damping of the least damped mode present in the data', fontsize=FS_TIT)
ax.legend(fontsize=FS_LEG)
ax.grid(True, which='both', alpha=0.3)
fig.savefig('predicted_rho01.png', bbox_inches='tight')
plt.close(fig)


# ---------------------------------------------------------------------------
# Iteration counts vs dt / vs the groups.
# ---------------------------------------------------------------------------
def counts_plot(xfun, xlabel, fname, ykey, ylabel, vline=None):
    fig, ax = plt.subplots(figsize=(10, 6))
    any_data = False
    for mem, per_dt in zip(members, records):
        x = [xfun(per_dt[d]) for d in dts]
        y = [per_dt[d][ykey] for d in dts]
        if np.all(np.isnan(y)):
            continue
        any_data = True
        ax.semilogx(x, y, marker='o', label=mem['legend'])
    if not any_data:
        plt.close(fig)
        return
    if vline is not None:
        ax.axvline(vline, ls=':', color='grey')
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_xlabel(xlabel, fontsize=FS_LAB)
    ax.set_ylabel(ylabel, fontsize=FS_LAB)
    ax.set_title(f'{ylabel} vs {xlabel} for {SWEEP_TITLE}', fontsize=FS_TIT)
    ax.legend(fontsize=FS_LEG)
    ax.grid(True, which='both', alpha=0.3)
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)


xo = cps.crossover_dt(args.length, args.height, N) if test != 'ar' else None
if mode == 'linear':
    counts_plot(lambda r: r['dt'], r'$\Delta t$', f'error_{T_NAME}.png', 'its', 'KSP iterations', vline=xo)
    counts_plot(lambda r: r['group'], r'$\varepsilon\,\Delta t\,N$', f'error_{T_NAME}_vs_group.png', 'its',
                'KSP iterations', vline=1.0)
    counts_plot(lambda r: r['g4'], r'$4\varepsilon^2 m$', f'error_{T_NAME}_vs_4eps2m.png', 'its',
                'KSP iterations', vline=1.0)
else:
    counts_plot(lambda r: r['dt'], r'$\Delta t$', f'snes_ksp_{T_NAME}.png', 'total_ksp',
                'total KSP iterations', vline=xo)
    counts_plot(lambda r: r['group'], r'$\varepsilon\,\Delta t\,N$', f'snes_ksp_{T_NAME}_vs_group.png', 'total_ksp',
                'total KSP iterations', vline=1.0)
    counts_plot(lambda r: r['g4'], r'$4\varepsilon^2 m$', f'snes_ksp_{T_NAME}_vs_4eps2m.png', 'total_ksp',
                'total KSP iterations', vline=1.0)
    counts_plot(lambda r: r['dt'], r'$\Delta t$', f'newton_steps_{T_NAME}.png', 'newton', 'Newton steps', vline=xo)
    if test == 'dt' and os.path.exists('snes_ksp_shift.png'):
        shutil.copyfile('snes_ksp_shift.png', 'error_shift_nonlinear.png')


# ---------------------------------------------------------------------------
# dt test: the (dt, K) map with the analytical lines K = 1 and K_operative(dt).
# ---------------------------------------------------------------------------
def _log_edges(x):
    x = np.asarray(sorted(x), dtype=float)
    if len(x) == 1:
        return np.array([x[0] / np.sqrt(10), x[0] * np.sqrt(10)])
    mids = np.sqrt(x[:-1] * x[1:])
    return np.concatenate([[x[0] ** 2 / mids[0]], mids, [x[-1] ** 2 / mids[-1]]])


if test == 'dt':
    Ks = [m['K'] for m in members]
    Z = np.full((len(Ks), len(dts)), np.nan)
    for i, per_dt in enumerate(records):
        for j, d in enumerate(dts):
            Z[i, j] = per_dt[d][count_key]
    if np.any(np.isfinite(Z)):
        fig, ax = plt.subplots(figsize=(10, 6))
        pc = ax.pcolormesh(_log_edges(dts), _log_edges(Ks), np.ma.masked_invalid(Z), cmap='viridis', shading='flat')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ! hatch the cells that did not converge / diverged / are missing
        for i, per_dt in enumerate(records):
            for j, d in enumerate(dts):
                r = per_dt[d]
                if r['status'] != 'ok':
                    ex, ey = _log_edges(dts), _log_edges(Ks)
                    ax.add_patch(matplotlib.patches.Rectangle((ex[j], ey[i]), ex[j + 1] - ex[j], ey[i + 1] - ey[i],
                                                              fill=False, hatch='//', edgecolor='red', lw=0))
                    ax.text(np.sqrt(ex[j] * ex[j + 1]), np.sqrt(ey[i] * ey[i + 1]), CELL.get(r['status'], r['status']),
                            ha='center', va='center', fontsize=7, color='red')
                elif np.isfinite(r[count_key]):
                    ex, ey = _log_edges(dts), _log_edges(Ks)
                    ax.text(np.sqrt(ex[j] * ex[j + 1]), np.sqrt(ey[i] * ey[i + 1]), f"{int(r[count_key])}",
                            ha='center', va='center', fontsize=8, color='white')
        ax.plot(dt_fine, [cps.K_operative(d, args.length, args.height, N) for d in dt_fine], 'w--', lw=2,
                label=r'$K_{op} = 1/\min(1, 4\varepsilon^2 m)$ ($\delta = \sigma_{01}/2$)')
        ax.axhline(1.0, color='w', ls=':', lw=2, label=r'$K = 1$ ($\delta = \delta_{thr}$)')
        ax.axvline(cps.crossover_dt(args.length, args.height, N), color='grey', ls=':', lw=1)
        ax.set_xlim(_log_edges(dts)[[0, -1]])
        ax.set_ylim(_log_edges(Ks)[[0, -1]])
        ax.set_xlabel(r'$\Delta t$ (s)', fontsize=FS_LAB)
        ax.set_ylabel(r'$K = \delta / \delta_{thr}$', fontsize=FS_LAB)
        ax.set_title(('KSP iterations' if mode == 'linear' else 'total KSP iterations') + r' over the $(\Delta t, K)$ plane',
                     fontsize=FS_TIT)
        fig.colorbar(pc, ax=ax, label='iterations')
        ax.legend(fontsize=11, loc='upper left')
        fig.savefig(f'iteration_map_{T_NAME}.png', bbox_inches='tight')
        plt.close(fig)


# ---------------------------------------------------------------------------
# Convergence histories at fixed dt, one curve per swept value.
# ---------------------------------------------------------------------------
def history_at_fixed_dt(dt):
    if mode == 'linear':
        fig, ax = plt.subplots(figsize=(10, 6))
        any_data = False
        for mem, per_dt in zip(members, records):
            err = per_dt[dt]['error']
            if err is None or len(err) == 0:
                continue
            ax.semilogy(np.arange(len(err)), np.clip(err, 1e-16, None), marker='o', label=mem['legend'])
            any_data = True
        if any_data:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            ax.set_xlabel('KSP iteration', fontsize=FS_LAB)
            ax.set_ylabel('relative error', fontsize=FS_LAB)
            ax.set_title(rf'Error convergence at $\Delta t$={dt:g} s, {SWEEP_TITLE}', fontsize=FS_TIT)
            ax.legend(fontsize=FS_LEG)
            ax.grid(True, which='both', alpha=0.3)
            fig.savefig(f'error_{T_NAME}_Robust_dt{dt}.png', bbox_inches='tight')
        plt.close(fig)
    else:
        fig, ax = plt.subplots(figsize=(10, 6))
        any_data = False
        for mem, per_dt in zip(members, records):
            se, cum = per_dt[dt]['snes_error'], per_dt[dt]['snes_cum']
            if se is None or cum is None or len(se) == 0 or len(se) != len(cum):
                continue
            if len(se) > 1:   # ! last point is the converged state, as in plot_shift.py
                se, cum = se[:-1], cum[:-1]
            ax.semilogy(cum, np.clip(se, 1e-16, None), marker='o', label=mem['legend'])
            any_data = True
        if any_data:
            ax.set_xlabel('cumulative KSP iterations', fontsize=FS_LAB)
            ax.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$', fontsize=FS_LAB)
            ax.set_title(rf'SNES error vs cumulative KSP iterations at $\Delta t$={dt:g} s, {SWEEP_TITLE}',
                         fontsize=FS_TIT)
            ax.legend(fontsize=FS_LEG)
            ax.grid(True, which='both', alpha=0.3)
            fig.savefig(f'snes_error_{T_NAME}_dt{dt}.png', bbox_inches='tight')
        plt.close(fig)


for dt in dts:
    history_at_fixed_dt(dt)

# ! dt test, linear: the complementary view, history at fixed K over all dt.
if test == 'dt' and mode == 'linear':
    for mem, per_dt in zip(members, records):
        fig, ax = plt.subplots(figsize=(10, 6))
        any_data = False
        for dt in dts:
            err = per_dt[dt]['error']
            if err is None or len(err) == 0:
                continue
            ax.semilogy(np.arange(len(err)), np.clip(err, 1e-16, None), marker='o', label=rf'$\Delta t$={_sci(dt)}')
            any_data = True
        if any_data:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
            ax.set_xlabel('KSP iteration', fontsize=FS_LAB)
            ax.set_ylabel('relative error', fontsize=FS_LAB)
            ax.set_title(f"Error convergence at {mem['legend']}", fontsize=FS_TIT)
            ax.legend(fontsize=FS_LEG)
            ax.grid(True, which='both', alpha=0.3)
            fig.savefig(f"error_Robust_K{mem['K']!r}.png", bbox_inches='tight')
        plt.close(fig)

counts = {}
for per_dt in records:
    for d in dts:
        counts[per_dt[d]['status']] = counts.get(per_dt[d]['status'], 0) + 1
print(f"{mode} {test}: " + ', '.join(f'{v} {k}' for k, v in sorted(counts.items())) + f", out of {len(members) * len(dts)} runs.")
print(open('summary_table.md').read())
