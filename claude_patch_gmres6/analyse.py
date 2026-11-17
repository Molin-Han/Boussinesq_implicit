'''
Figures, tables and per-folder READMEs of the 2026-11-16 study: star patch vs single-column (vanka) patch,
both codes, GMRES(6) smoothing, the paper's robustness sweeps (delta = C1/dt), plus an MPI-robustness study.

    python analyse.py [--only main paper smoother repeat mpi overview]

Inputs (read only):
    ../../data/patch_gmres6/data_<code>_<test>_deg2[_ew]_gmres6_<patch>/   the new sweeps (2 MPI ranks)
    ../../data/patch_gmres6/mpi/<code>_<patch>/np<NN>_r<R>/ + points.csv   the MPI study (make_mpi_jobs.py)
    ../../data/paper_robustness/data_<code>_<test>_deg2[_ew]_<patch>/      the 2026-10-16 runs (linear: Chebyshev(6),
                                                                            outer GMRES-left; nonlinear: GMRES(6), identical
                                                                            configuration to the new runs)
Outputs: figures/<folder>/ (PNG, CSV, README.md); see ../claude_patch_gmres6/README.md.

Counting conventions (the same as claude_small_patch/compare_paper_patches.py):
    linear     n = number of outer FGMRES iterations of the monitored solve = lines of error_*.out - 1
    nonlinear  total KSP (FGMRES) iterations over all Newton steps of the monitored time step = last entry of
               snes_ksp_cum_*.out; Newton steps = entries - 1
    status     ok | maxit (a linear solve, or one Newton step's KSP, reached --maxit 100) | div (no output file and the
               job failed: ConvergenceError / DIVERGED_*) | crash / timeout | missing (not run yet)
'''
import argparse
import csv
import glob
import math
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
YR2 = os.path.abspath(os.path.join(REPO, '..'))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, 'claude_small_patch'))
import claude_paper_shift as cps           # noqa: E402
import paper_style as ps                   # noqa: E402
import restyle_paper_figures as R          # noqa: E402  (plotting functions only; nothing runs on import)

NEW = os.path.join(YR2, 'data', 'patch_gmres6')
OLD = os.path.join(YR2, 'data', 'paper_robustness')
MPI = os.path.join(NEW, 'mpi')
FIG = os.path.join(HERE, 'figures')
MAXIT = 100
LX, N_FREQ = 3.0e5, 1.0e-2
TESTS = ('dt', 'ar', 'dx', 'dz')
CODES = ('linear', 'nonlinear')
PATCHES = ('star', 'vanka')
PNAME = {'star': 'star', 'vanka': 'single column'}
SWEEP = {
    'linear': dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000.], C1=1e-3,
                   C1_list=[1., 0.1, 0.01, 1e-3, 1e-4, 1e-5], heights=[12000, 8000, 4000, 2000, 1000, 500, 250],
                   nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80]),
    'nonlinear': dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000., 50000.], C1=1e-5,
                      C1_list=[1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6], heights=[12000, 8000, 4000, 2000, 1000],
                      nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80]),
}
HIST_DT = {'linear': {'dt': 50., 'ar': 50., 'dx': 50., 'dz': 50.},          # the dt of the paper's history figures
           'nonlinear': {'dt': 500., 'ar': 100., 'dx': 50., 'dz': 50.}}
CONFIG = {
    'linear': ('LB_time_slice.py --richardson --reordering --rtol 1e-8 --maxit 100 --mg_smoother gmres --patch <p>: '
               'outer FGMRES (right, true residual, rtol 1e-8); Schur fieldsplit with exact pressure-mass solve; '
               'shifted auxiliary operator for the velocity-buoyancy block solved by ONE full-multigrid cycle '
               '(3 levels, horizontal coarsening, MUMPS/LU coarse); level smoother GMRES(6) + additive Schwarz on '
               'column patches (ASMStarPC / ASMVankaPC, LU with RCM); f = 0; 2 MPI ranks.'),
    'nonlinear': ('Nonlinear_Boussinesq_Irk_SC.py --no_rotation -firedrake_1_snes_ksp_ew -firedrake_1_snes_ksp_ew_rtol0 1e-2 '
                  '-firedrake_1_snes_atol 1e-09 --maxit 100 --patch <p> (its default solver): Newton (newtonls) with '
                  'Eisenstat-Walker (PETSc v2, rtol_0 = 1e-2): the FGMRES tolerance of every Newton step is set by EW, so the '
                  "script's --rtol 1e-6 is overridden; outer FGMRES (right, atol 1e-9, at most 100 iterations per Newton "
                  'step, restart 30); IRKAuxiliaryOperatorPC with the shifted '
                  'operator (pressure elimination p -> p - div(u)/delta); Schur fieldsplit, exact pressure-mass solve; ONE '
                  'full-multigrid cycle on the velocity-buoyancy block with GMRES(6) + ASM column-patch smoothing, '
                  'MUMPS/LU coarse; 2 MPI ranks.'),
}


# --------------------------------------------------------------------------------------------- data access
def folder_new(code, test, patch):
    return os.path.join(NEW, f'data_{code}_{test}_deg2' + ('_ew' if code == 'nonlinear' else '') + f'_gmres6_{patch}')


def folder_old(code, test, patch):
    return os.path.join(OLD, f'data_{code}_{test}_deg2' + ('_ew' if code == 'nonlinear' else '') + f'_{patch}')


def members(code, test):
    S = SWEEP[code]
    return {'dt': S['C1_list'], 'ar': S['heights'], 'dx': S['nxs'], 'dz': S['nzs']}[test]


def key_of(test):
    return {'dt': 'shift', 'ar': 'ar', 'dx': 'nx', 'dz': 'nz'}[test]


def target(test, val, dt):
    return {'dt': val / dt, 'ar': val / LX}.get(test, val)


def geometry(code, test, val, dt):
    '''(nx, nz, height, shift) of a sweep point.'''
    nx, nz, h = 20, 20, 4000.
    if test == 'ar':
        h = float(val)
    elif test == 'dx':
        nx = int(val)
    elif test == 'dz':
        nz = int(val)
    shift = val / dt if test == 'dt' else SWEEP[code]['C1'] / dt
    return nx, nz, h, shift


def groups(code, test, val, dt):
    nx, nz, h, shift = geometry(code, test, val, dt)
    K = shift / cps.shift_threshold(dt, LX, h, N_FREQ)
    rho1 = cps.m_factor(dt, N_FREQ) * h ** 2 / (math.pi ** 2 * (LX / (4 * nx)) ** 2)
    return dict(nx=nx, nz=nz, height=h, shift=shift, K=K, rho1=rho1,
                eps_dt_N=cps.group_eps_dt_N(dt, LX, h, N_FREQ), regime=cps.regime(dt, LX, h, N_FREQ))


def find(d, prefix, dt, key, tgt):
    best, bd = None, 1e9
    for f in glob.glob(f'{d}/{prefix}_dt{float(dt)!r}_{key}*.out'):
        m = re.search(rf'_{key}(.+)\.out$', os.path.basename(f))
        try:
            v = float(m.group(1))
        except (AttributeError, ValueError):
            continue
        dd = abs(np.log(v) - np.log(tgt)) if (v > 0 and tgt > 0) else abs(v - tgt)
        if dd < bd:
            best, bd = f, dd
    return best if bd <= 0.2 else None


def job_table(d):
    '''status of every job of a data folder, keyed by its parameters.'''
    out = []
    for st in glob.glob(f'{d}/jobs/*/status'):
        try:
            rc, secs, reason = open(st).read().split()[:3]
            cmd = open(os.path.join(os.path.dirname(st), 'cmd.sh')).read()
        except (OSError, ValueError):
            continue
        g = lambda k, cast=float: cast(re.search(rf'--{k} (\S+)', cmd).group(1))   # noqa: E731
        out.append(dict(dt=g('dt'), shift=g('shift'), height=g('height'), nx=g('nx', int), nz=g('nz', int),
                        rc=int(rc), secs=int(secs), reason=reason, dir=os.path.dirname(st)))
    return out


def final_rel_residual(jobdir):
    '''last ||r(i)||/||b|| printed by -ksp_monitor_true_residual in the job's run.log (the monitored solve).'''
    try:
        vals = re.findall(r'\|\|r\(i\)\|\|/\|\|b\|\| ([0-9.eE+-]+)', open(os.path.join(jobdir, 'run.log')).read())
        return float(vals[-1]) if vals else np.nan
    except OSError:
        return np.nan


LINEAR_RTOL = 1e-8


def match_job(jobs, code, test, val, dt):
    nx, nz, h, shift = geometry(code, test, val, dt)
    for j in jobs:
        if abs(j['dt'] - dt) <= 1e-9 * max(1., dt) and j['nx'] == nx and j['nz'] == nz and abs(j['height'] - h) < 1e-6 \
                and abs(math.log(j['shift']) - math.log(shift)) < 0.2:
            return j
    return None


def record(d, code, test, val, dt, jobs=None):
    rec = dict(count=np.nan, newton=np.nan, status='missing', error=None, cum=None, secs=np.nan)
    key, tgt = key_of(test), target(test, val, dt)
    err = find(d, 'error', dt, key, tgt)
    if err is not None:
        rec['error'] = np.atleast_1d(np.loadtxt(err))
    if code == 'linear':
        if err is not None:
            n = len(rec['error']) - 1
            rec['status'] = 'maxit' if n >= MAXIT else 'ok'
            rec['count'] = n
    else:
        cf = find(d, 'snes_ksp_cum', dt, key, tgt)
        if cf is not None:
            cum = np.atleast_1d(np.loadtxt(cf)).astype(int)
            rec['cum'] = cum
            rec['newton'] = len(cum) - 1
            per = np.diff(cum) if len(cum) > 1 else cum
            rec['status'] = 'maxit' if np.any(per >= MAXIT) else 'ok'
            rec['count'] = int(cum[-1])
            se = find(d, 'snes_error', dt, key, tgt)
            rec['snes_error'] = np.atleast_1d(np.loadtxt(se)) if se else None
    jobs = job_table(d) if jobs is None else jobs
    j = match_job(jobs, code, test, val, dt)
    if j is not None:
        rec['secs'] = j['secs']
        # ! with ksp_converged_maxits PETSc reports CONVERGED_ITS at the cap even when rtol is met at exactly
        # ! iteration 100 (it tests max_it first): such a linear run converged, it did not fail.
        if code == 'linear' and rec['status'] == 'maxit' and final_rel_residual(j['dir']) <= LINEAR_RTOL:
            rec['status'] = 'ok'
            rec['at_cap'] = True
        if rec['status'] == 'missing' and j['rc'] != 0:
            r = j['reason']
            rec['status'] = 'div' if ('Converg' in r or 'DIVERGED' in r) else ('timeout' if r == 'timeout' else 'crash')
    return rec


def ok(r):
    return r['status'] == 'ok'


def ratio(a, b):
    return a['count'] / b['count'] if ok(a) and ok(b) and b['count'] > 0 else np.nan


def cell(r):
    if ok(r):
        return f"{int(r['count'])}" + (f"/{int(r['newton'])}" if np.isfinite(r['newton']) else '')
    return r['status'] if r['status'] != 'maxit' else f">={MAXIT}"


def member_label(code, test, val):
    return {'dt': f'C1={val:g}', 'ar': f'Lz={val / 1000:g} km', 'dx': f'Nx={val}', 'dz': f'Nz={val}'}[test]


def member_tex(code, test, val):
    return {'dt': rf'$C_1 = {ps.sci(val)[1:-1]}$', 'ar': rf'$L_z = {val / 1000:g}$ km',
            'dx': rf'$N_x = {val}$', 'dz': rf'$N_z = {val}$'}[test]


def collect(code, test):
    '''all records of one sweep: rows = members x dts.'''
    jn = {p: job_table(folder_new(code, test, p)) for p in PATCHES}
    jo = {p: job_table(folder_old(code, test, p)) for p in PATCHES}
    rows = []
    for val in members(code, test):
        for dt in SWEEP[code]['dts']:
            row = dict(code=code, test=test, member=val, dt=dt, **groups(code, test, val, dt))
            for p in PATCHES:
                row[p] = record(folder_new(code, test, p), code, test, val, dt, jn[p])
                row['old_' + p] = record(folder_old(code, test, p), code, test, val, dt, jo[p])
            rows.append(row)
    return rows


# --------------------------------------------------------------------------------------------- helpers
def write_csv(path, rows, cols):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, '') for c in cols])


def fmt(x, nd=3):
    if isinstance(x, (float, np.floating)):
        if not np.isfinite(x):
            return ''
        return f'{x:.{nd}g}'
    return str(x)


def stats_line(vals):
    v = np.array([x for x in vals if np.isfinite(x)])
    if not len(v):
        return 'no data'
    return f'n = {len(v)}, median {np.median(v):.2f}, mean {np.mean(v):.2f}, min {v.min():.2f}, max {v.max():.2f}'


def status_counts(rows, k):
    c = {}
    for r in rows:
        c[r[k]['status']] = c.get(r[k]['status'], 0) + 1
    return ', '.join(f'{s}: {n}' for s, n in sorted(c.items()))


def md_table(header, body):
    out = ['| ' + ' | '.join(header) + ' |', '|' + '---|' * len(header)]
    out += ['| ' + ' | '.join(str(x) for x in row) + ' |' for row in body]
    return '\n'.join(out)


# --------------------------------------------------------------------------------------------- 1. main comparison folders
def main_folder(code, test):
    rows = collect(code, test)
    name = f'{"1" if code == "linear" else "2"}{TESTS.index(test) + 1}_{code}_{test}'
    out = os.path.join(FIG, name)
    os.makedirs(out, exist_ok=True)
    mem = members(code, test)
    dts = SWEEP[code]['dts']
    colors = ps.ramp(len(mem))

    # counts against dt, star | single column | ratio
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.8))
    for i, val in enumerate(mem):
        for ax, p in zip(axs[:2], PATCHES):
            y = [r[p]['count'] if ok(r[p]) else np.nan for r in rows if r['member'] == val]
            ax.plot(dts, y, color=colors[i], **(ps.STAR if p == 'star' else ps.COLUMN), label=member_tex(code, test, val))
            for r in rows:
                if r['member'] == val and r[p]['status'] in ('maxit', 'div', 'crash', 'timeout'):
                    ax.plot(r['dt'], MAXIT * 1.08, marker='x', color=colors[i], ms=6, clip_on=False)
        y = [ratio(r['vanka'], r['star']) for r in rows if r['member'] == val]
        axs[2].plot(dts, y, color=colors[i], marker='o', label=member_tex(code, test, val))
    ylab = 'outer FGMRES iterations' if code == 'linear' else 'total FGMRES iterations (all Newton steps)'
    for ax, p in zip(axs[:2], PATCHES):
        ax.set_xscale('log')
        ax.set_ylim(0, MAXIT * 1.12)
        ps.finish(ax, r'$\Delta t$ (s)', ylab, legend=(p == 'star'), small=True, title=f'{PNAME[p]} patch')
    axs[2].set_xscale('log')
    axs[2].axhline(1.0, color=ps.GREY, lw=0.8, ls=':')
    ps.finish(axs[2], r'$\Delta t$ (s)', 'single column / star', legend=False, title='iteration ratio')
    fig.suptitle(f'{code} code, {test} sweep, GMRES(6) smoothing (x at the top: maxit / diverged)', fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'counts_vs_dt.png'), dpi=160)
    plt.close(fig)

    # ratio map with the counts written in the cells
    R_ = np.array([[ratio(r['vanka'], r['star']) for r in rows if r['member'] == val] for val in mem])
    fig, ax = plt.subplots(figsize=(1.0 + 0.9 * len(dts), 0.9 + 0.45 * len(mem)))
    from matplotlib.colors import TwoSlopeNorm
    im = ax.imshow(np.log2(R_), cmap='RdBu_r', norm=TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=3.5), aspect='auto')
    for i, val in enumerate(mem):
        for k, dt in enumerate(dts):
            r = [x for x in rows if x['member'] == val and x['dt'] == dt][0]
            ax.text(k, i, f"{cell(r['star']).split('/')[0]}\n{cell(r['vanka']).split('/')[0]}", ha='center', va='center', fontsize=6.5)
    ax.set_xticks(range(len(dts)), [f'{d:g}' for d in dts])
    ax.set_yticks(range(len(mem)), [member_label(code, test, v) for v in mem])
    ax.set_xlabel(r'$\Delta t$ (s)')
    cb = fig.colorbar(im, ax=ax, fraction=0.04)
    cb.set_label('log2(single column / star)')
    ax.set_title('cell: star count (top) / single-column count (bottom)', fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'ratio_map.png'), dpi=160)
    plt.close(fig)

    # error histories, star (solid) vs single column (dashed)
    # linear: ||x_k - x*||/||x*|| per FGMRES iteration (error_*.out); nonlinear: error of each Newton iterate
    # (snes_error_*.out) against the cumulative FGMRES count at that iterate (snes_ksp_cum_*.out), as in the
    # paper's nonlinear figures.  (The nonlinear error_*.out is not used: it compares every FGMRES iterate with the
    # correction of the warm-up pass's LAST Newton step, so only its last segment is a convergence history.)
    hdt = HIST_DT[code][test]
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    for i, val in enumerate(mem):
        r = [x for x in rows if x['member'] == val and x['dt'] == hdt][0]
        for p in PATCHES:
            if code == 'linear':
                e = r[p]['error']
                xs = np.arange(len(e)) if e is not None else None
            else:
                e, xs = r[p].get('snes_error'), r[p].get('cum')
            if e is not None and xs is not None and len(e) > 1 and len(xs) == len(e):
                ax.semilogy(xs, np.maximum(e, 1e-16), color=colors[i], lw=1.3, marker='o' if code == 'nonlinear' else None,
                            ms=3, ls='-' if p == 'star' else '--', label=member_tex(code, test, val) if p == 'star' else None)
    ps.finish(ax, 'cumulative FGMRES iterations (markers: Newton iterates)' if code == 'nonlinear' else 'FGMRES iteration',
              r'$\|u_k - u^*\| / \|u^*\|$ (Newton iterates)' if code == 'nonlinear' else r'$\|x_k - x^*\| / \|x^*\|$',
              small=True, title=rf'$\Delta t = {hdt:g}$ s; solid star, dashed single column')
    fig.tight_layout()
    fig.savefig(os.path.join(out, f'histories_dt{hdt:g}.png'), dpi=160)
    plt.close(fig)

    # csv
    flat = []
    for r in rows:
        d = {k: r[k] for k in ('code', 'test', 'member', 'dt', 'nx', 'nz', 'height', 'shift', 'K', 'rho1', 'eps_dt_N', 'regime')}
        for p in PATCHES:
            d[f'{p}_count'] = fmt(r[p]['count'])
            d[f'{p}_newton'] = fmt(r[p]['newton'])
            d[f'{p}_status'] = r[p]['status']
            d[f'{p}_secs'] = fmt(r[p]['secs'])
            d[f'{p}_converged_at_cap'] = 'yes' if r[p].get('at_cap') else ''
            d[f'old_{p}_count'] = fmt(r['old_' + p]['count'])
            d[f'old_{p}_status'] = r['old_' + p]['status']
        d['ratio_vanka_star'] = fmt(ratio(r['vanka'], r['star']))
        flat.append(d)
    cols = list(flat[0].keys())
    write_csv(os.path.join(out, 'counts.csv'), flat, cols)

    # README
    rat = [ratio(r['vanka'], r['star']) for r in rows]
    fails_s = {(r['member'], r['dt']) for r in rows if r['star']['status'] not in ('ok', 'missing')}
    fails_v = {(r['member'], r['dt']) for r in rows if r['vanka']['status'] not in ('ok', 'missing')}
    worst = sorted([r for r in rows if np.isfinite(ratio(r['vanka'], r['star']))], key=lambda r: -ratio(r['vanka'], r['star']))[:8]
    best = sorted([r for r in rows if np.isfinite(ratio(r['vanka'], r['star']))], key=lambda r: ratio(r['vanka'], r['star']))[:5]
    bands = []
    for lo, hi, lab in ((0, 0.1, 'K < 0.1'), (0.1, 1, '0.1 <= K < 1'), (1, 10, '1 <= K < 10'), (10, np.inf, 'K >= 10')):
        bands.append((lab, stats_line([ratio(r['vanka'], r['star']) for r in rows if lo <= r['K'] < hi])))
    table = []
    for val in mem:
        table.append([member_label(code, test, val)] + [f"{cell(r['star'])} / {cell(r['vanka'])}" for r in rows if r['member'] == val])
    old_note = ('Chebyshev(6) smoothing with outer GMRES-left (preconditioned residual), so they are NOT the same solver; '
                'see ../30_smoother_linear/' if code == 'linear' else
                'the same solver configuration (GMRES(6), FGMRES-right, EW) run on 2026-10-16; differences are run-to-run '
                'noise, see ../40_repeatability_nonlinear/')
    txt = f'''# {code} code, {test} sweep: star vs single-column patch, GMRES(6) smoothing

Generated by `claude_patch_gmres6/analyse.py`; do not edit by hand (re-run the script).

## Configuration
{CONFIG[code]}

Sweep: the paper's {test} sweep (`Bash_Generator/bash_paper_patch.py --mode {code} --test {test}`), delta = C1/dt with
{'C1 in ' + str(SWEEP[code]['C1_list']) if test == 'dt' else 'C1 = ' + str(SWEEP[code]['C1'])}; members {member_label(code, test, mem[0])} ... {member_label(code, test, mem[-1])};
dt in {dts} s.

Data (all paths relative to `yr2/`): `{os.path.relpath(folder_new(code, test, 'star'), YR2)}` and
`{os.path.relpath(folder_new(code, test, 'vanka'), YR2)}` (each run's command in `jobs/job_*/cmd.sh`, its log in
`jobs/job_*/run.log`, status in `jobs/job_*/status`; take C1 from the `# C1=` comment of cmd.sh or from counts.csv, not
from run.log, which prints shift*dt^1.5 under that name). Columns `old_*` in counts.csv are the 2026-10-16 runs in
`{os.path.relpath(folder_old(code, test, 'star'), YR2)}` etc.: {old_note}.

Counting: {'n = outer FGMRES iterations to rtol 1e-8 on the true residual; a run that meets rtol at exactly iteration 100 counts as converged (column *_converged_at_cap)' if code == 'linear' else 'total FGMRES iterations over all Newton steps of the monitored time step (counts.csv also gives Newton steps)'};
`>=100` = a solve hit maxit 100 without converging; `div` = Newton/KSP diverged (job failed without output).
K = delta / delta_thr(dt) (claude_paper_shift.py): K < 1 is the regime where the grad-div term dt/delta dominates; it
occurs at large dt and small C1 (in these sweeps K does not depend on Lz for dt <= 1000 s, the weakly stratified branch).
rho1 = m Lz^2 / (pi^2 dx_fine^2), m = 1 + (dt N)^2/4, dx_fine = Lx / (4 Nx) (the 2026-09-29 study's kernel parameter).

## Figures
- `counts_vs_dt.png`: iterations against dt, one curve per member, star | single column | ratio (single column / star).
- `ratio_map.png`: log2 ratio per (member, dt), colour range -1 ... 3.5 (ratio 0.5 ... 11); each cell gives the two counts (star on top).
- `histories_dt{hdt:g}.png`: at dt = {hdt:g} s{'' if code == 'nonlinear' or test != 'dt' else ' (the paper shows the shift sweep as counts against dt, not as histories)'},
  star solid, single column dashed. ''' + ('''Relative error ||x_k - x*||/||x*|| against FGMRES iteration, where x* is the
  result of an identical warm-up solve in the same run (same iterative solver and tolerance), not a direct solve: the final
  drop to about 1e-15 is an artefact, the penultimate value measures the last update.''' if code == 'linear' else '''Relative
  error of each Newton iterate (snes_error_*.out, against the warm-up pass's converged state) against the cumulative FGMRES
  count at that iterate (snes_ksp_cum_*.out), one marker per Newton iterate, as in the paper's nonlinear figures.''') + f'''
- `counts.csv`: every point with K, rho1, eps*dt*N, both patches' counts, Newton steps, status, run time, and the old runs.

## Numbers (auto-generated)
- runs: star [{status_counts(rows, 'star')}]; single column [{status_counts(rows, 'vanka')}]
- failure sets (maxit/div/crash) identical: {fails_s == fails_v}; star only: {sorted(fails_s - fails_v)}; single column only: {sorted(fails_v - fails_s)}
- ratio single column / star over all points where both converged: {stats_line(rat)}
- by K: ''' + '; '.join(f'{lab}: {s}' for lab, s in bands) + f'''
- largest ratios: ''' + '; '.join(f"{member_label(code, test, r['member'])}, dt={r['dt']:g}: {cell(r['star'])} vs {cell(r['vanka'])} (K={r['K']:.2g})" for r in worst) + f'''
- smallest ratios: ''' + '; '.join(f"{member_label(code, test, r['member'])}, dt={r['dt']:g}: {cell(r['star'])} vs {cell(r['vanka'])}" for r in best) + f'''

## Table: star / single column (count{'/Newton steps' if code == 'nonlinear' else ''}) per member (rows) and dt (columns)
{md_table(['member'] + [f'dt={d:g}' for d in dts], table)}
'''
    open(os.path.join(out, 'README.md'), 'w').write(txt)
    print('wrote', out)
    return rows


# --------------------------------------------------------------------------------------------- 2. paper-style figures
def paper_folder():
    ps.apply()
    for p, sub in (('star', 'star'), ('vanka', 'single_column')):
        for corrected in (False, True):
            out = os.path.join(FIG, '50_paper_style', sub, 'corrected' if corrected else 'paper_conventions')
            os.makedirs(out, exist_ok=True)
            f = lambda c, t: folder_new(c, t, p)   # noqa: E731
            jobs = [
                (R.linear_history, f('linear', 'dx'), 'dx', 50.0, corrected, out, 'error_dx_Robust_dt50.0.png', ps.SINGLE),
                (R.linear_history, f('linear', 'dz'), 'dz', 50.0, corrected, out, 'error_dz_Robust_dt50.0.png', ps.SINGLE),
                (R.linear_history, f('linear', 'ar'), 'ar', 50.0, corrected, out, 'error_AR_Robust_dt50.0.png', ps.SINGLE),
                (R.linear_shift_sweep, f('linear', 'dt'), corrected, out, 'error_shift.png', ps.SINGLE),
                (R.nonlinear_history, f('nonlinear', 'dx'), 'dx', 50.0, corrected, out, 'snes_error_dx_dt50.0.png', ps.HALF),
                (R.nonlinear_history, f('nonlinear', 'dz'), 'dz', 50.0, corrected, out, 'snes_error_dz_dt50.0.png', ps.HALF),
                (R.nonlinear_history, f('nonlinear', 'ar'), 'ar', 100.0, corrected, out, 'snes_error_AR_dt100.0.png', ps.HALF),
                (R.nonlinear_history, f('nonlinear', 'ar'), 'ar', 5000.0, corrected, out, 'snes_error_AR_dt5000.0.png', ps.HALF, True),
                (R.nonlinear_history, f('nonlinear', 'dt'), 'dt', 500.0, corrected, out, 'snes_error_shift_dt500.0.png', ps.HALF),
                (R.nonlinear_shift_sweep, f('nonlinear', 'dt'), corrected, out, 'error_shift_nonlinear.png', ps.HALF),
            ]
            for fn, *args in jobs:
                try:
                    fn(*args)
                except Exception as e:          # partial data while the sweeps are running
                    print('  skipped', args[-2] if isinstance(args[-2], str) else '', type(e).__name__, e)
            plt.close('all')
    readme = '''# The paper's ten iteration figures, redrawn from the GMRES(6) runs

Generated by `claude_patch_gmres6/analyse.py` with the plotting functions of `claude_small_patch/restyle_paper_figures.py`
(same style, file names and panels as the paper's figures 4-12), from the new data in `data/patch_gmres6/`:

- `star/`           star patch (ASMStarPC)          = the paper's patch
- `single_column/`  single-column patch (ASMVankaPC)

Each has two versions:
- `paper_conventions/`: the paper's plotting conventions (linear histories drop the converged point; the linear shift
  sweep plots n + 1; the nonlinear shift sweep plots total KSP + Newton steps; the nonlinear histories keep the converged
  Newton iterate except in snes_error_AR_dt5000.0, as in the paper), i.e. drop-in replacements;
- `corrected/`: the shift sweeps plot the iteration counts themselves (n; total KSP over the Newton steps) and the linear
  histories keep their final point, but the nonlinear histories DROP the converged Newton iterate (it is the reference
  state), so their last abscissa is the cumulative count before the final Newton step, not the total; snes_error_AR_dt5000.0
  is therefore the same in both versions. For nonlinear totals use paper_conventions/ or counts.csv of ../2x_*.

Solvers. These runs: GMRES(6) level smoothing for BOTH codes and outer FGMRES with right preconditioning (the solver the
paper's Section 5.1 text describes). The paper's current figures were made with other solvers (from the code at the time,
git history of Boussinesq_implicit): linear figures 4-7 with Chebyshev(6) smoothing and outer GMRES with LEFT
preconditioning (`LB_time_slice.py --richardson`); nonlinear figures 8-12 (data of 2026-05-20, commit 7771540) with
Chebyshev(6) smoothing and outer GMRES with right preconditioning and Eisenstat-Walker. GMRES(6) + FGMRES became the
nonlinear default in commit 82cd238 (2026-07-08). So both codes' figures here differ from the paper's; see
../30_smoother_linear/ for the linear comparison. All paths relative to `yr2/`.
'''
    open(os.path.join(FIG, '50_paper_style', 'README.md'), 'w').write(readme)
    print('wrote', os.path.join(FIG, '50_paper_style'))


# --------------------------------------------------------------------------------------------- 3. smoother comparison (linear)
def smoother_folder(allrows):
    out = os.path.join(FIG, '30_smoother_linear')
    os.makedirs(out, exist_ok=True)
    flat, lines = [], []
    fig, axs = plt.subplots(1, 4, figsize=(15, 3.8))
    for ax, test in zip(axs, TESTS):
        rows = allrows[('linear', test)]
        for p, mk, col in (('star', 'o', 'C0'), ('vanka', 's', 'C3')):
            x = [r['old_' + p]['count'] for r in rows if ok(r['old_' + p]) and ok(r[p])]
            y = [r[p]['count'] for r in rows if ok(r['old_' + p]) and ok(r[p])]
            if p == 'vanka':
                ax.scatter(x, y, s=22, marker=mk, facecolors='none', edgecolors=col, label=PNAME[p])
            else:
                ax.scatter(x, y, s=10, marker=mk, color=col, label=PNAME[p])
            # one-sided failures: plotted at the cap (x on the axis of the configuration that failed)
            for r in rows:
                o, n = ok(r['old_' + p]), ok(r[p])
                if o and not n:
                    ax.plot(r['old_' + p]['count'], MAXIT + 3, marker='x', color=col, ms=6, clip_on=False)
                elif n and not o:
                    ax.plot(MAXIT + 3, r[p]['count'], marker='x', color=col, ms=6, clip_on=False)
        lim = [0, MAXIT]
        ax.plot(lim, lim, color=ps.GREY, lw=0.8, ls=':')
        ax.set_xlim(0, MAXIT + 5); ax.set_ylim(0, MAXIT + 5)
        ps.finish(ax, 'Chebyshev(6), GMRES-left (2026-10-16)', 'GMRES(6), FGMRES-right (new)', small=True, title=f'linear {test} sweep')
        for r in rows:
            flat.append(dict(test=test, member=r['member'], dt=r['dt'], K=fmt(r['K']),
                             cheb_star=cell(r['old_star']), gmres_star=cell(r['star']),
                             cheb_vanka=cell(r['old_vanka']), gmres_vanka=cell(r['vanka']),
                             ratio_cheb=fmt(ratio(r['old_vanka'], r['old_star'])), ratio_gmres=fmt(ratio(r['vanka'], r['star']))))
        lines.append(f"| {test} | {stats_line([ratio(r['old_vanka'], r['old_star']) for r in rows])} | "
                     f"{stats_line([ratio(r['vanka'], r['star']) for r in rows])} | "
                     f"{stats_line([r['star']['count'] / r['old_star']['count'] for r in rows if ok(r['star']) and ok(r['old_star'])])} | "
                     f"{stats_line([r['vanka']['count'] / r['old_vanka']['count'] for r in rows if ok(r['vanka']) and ok(r['old_vanka'])])} |")
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'gmres_vs_chebyshev_counts.png'), dpi=160)
    plt.close(fig)
    # ratio distributions
    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    data, labels = [], []
    for test in TESTS:
        rows = allrows[('linear', test)]
        for key, lab in (('old', 'Cheb'), ('new', 'GMRES')):
            v = [ratio(r['old_vanka'], r['old_star']) if key == 'old' else ratio(r['vanka'], r['star']) for r in rows]
            data.append([x for x in v if np.isfinite(x)] or [np.nan])
            labels.append(f'{test}\n{lab}')
    ax.boxplot(data, tick_labels=labels, showfliers=True)
    ax.axhline(1.0, color=ps.GREY, lw=0.8, ls=':')
    ps.finish(ax, None, 'single column / star', legend=False, title='linear code: iteration ratio, Chebyshev(6) vs GMRES(6) smoothing')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'ratio_distribution.png'), dpi=160)
    plt.close(fig)
    write_csv(os.path.join(out, 'smoother_counts.csv'), flat, list(flat[0].keys()))
    txt = f'''# Linear code: GMRES(6) smoothing (new) vs Chebyshev(6) smoothing (2026-10-16 runs)

Generated by `claude_patch_gmres6/analyse.py`.

The 2026-10-16 linear star/single-column runs (`data/paper_robustness/data_linear_<test>_deg2_{{star,vanka}}`, also the data of
the paper's supplement SM2) used `LB_time_slice.py --richardson`: Chebyshev(6) level smoothing and outer GMRES with LEFT
preconditioning (stopping test on the preconditioned residual). The new runs use GMRES(6) smoothing and outer FGMRES with
RIGHT preconditioning (true residual), which is the solver the paper's Section 5.1 describes. Both: one full-multigrid
cycle, rtol 1e-8, maxit 100, 2 MPI ranks. The comparison changes the smoother, the outer Krylov method, the
preconditioning side and the residual norm at once, so their effects are not separated; compare trends, not single counts.

## Figures and tables
- `gmres_vs_chebyshev_counts.png`: new count against old count for every point, per sweep (filled dots star, open squares
  single column); an x on the top edge = converged with Chebyshev only, on the right edge = converged with GMRES(6) only.
- `ratio_distribution.png`: distribution of single column / star under each smoother, per sweep.
- `smoother_counts.csv`: every point: counts (or status) for the four combinations and both ratios.

## Numbers (auto-generated)
| sweep | single/star, Chebyshev | single/star, GMRES(6) | star: GMRES/Cheb | single column: GMRES/Cheb |
|---|---|---|---|---|
''' + '\n'.join(lines) + '\n'
    open(os.path.join(out, 'README.md'), 'w').write(txt)
    print('wrote', out)


# --------------------------------------------------------------------------------------------- 4. repeatability (nonlinear)
def repeat_folder(allrows):
    out = os.path.join(FIG, '40_repeatability_nonlinear')
    os.makedirs(out, exist_ok=True)
    flat, lines = [], []
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.8))
    for ax, p in zip(axs, PATCHES):
        diffs = []
        for test in TESTS:
            rows = allrows[('nonlinear', test)]
            x = [r['old_' + p]['count'] for r in rows if ok(r['old_' + p]) and ok(r[p])]
            y = [r[p]['count'] for r in rows if ok(r['old_' + p]) and ok(r[p])]
            ax.scatter(x, y, s=12, label=test)
            for r in rows:
                if ok(r['old_' + p]) and ok(r[p]):
                    diffs.append(r[p]['count'] - r['old_' + p]['count'])
                flat.append(dict(test=test, patch=p, member=r['member'], dt=r['dt'], K=fmt(r['K']),
                                 new=cell(r[p]), old=cell(r['old_' + p]),
                                 diff=fmt(r[p]['count'] - r['old_' + p]['count']) if ok(r[p]) and ok(r['old_' + p]) else ''))
        ax.plot([0, MAXIT], [0, MAXIT], color=ps.GREY, lw=0.8, ls=':')
        ax.set_xlim(0, MAXIT); ax.set_ylim(0, MAXIT)
        ps.finish(ax, '2026-10-16 run', '2026-11-16 run (same command)', small=True, title=f'{PNAME[p]} patch')
        d = np.array(diffs)
        allr = [r for t in TESTS for r in allrows[('nonlinear', t)]]
        one_failed = sum(1 for r in allr if ok(r[p]) != ok(r['old_' + p]))
        both_failed_diff = sum(1 for r in allr if not ok(r[p]) and not ok(r['old_' + p]) and r[p]['status'] != r['old_' + p]['status'])
        if len(d):
            lines.append(f'| {PNAME[p]} | {len(d)} | {int(np.sum(d == 0))} ({np.mean(d == 0) * 100:.0f}%) | {np.abs(d).max()} | '
                         f'{np.median(np.abs(d)):.1f} | {one_failed} | {both_failed_diff} |')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'new_vs_old_counts.png'), dpi=160)
    plt.close(fig)
    write_csv(os.path.join(out, 'repeat_counts.csv'), flat, list(flat[0].keys()))
    txt = '''# Nonlinear code: run-to-run repeatability (2026-11-16 vs 2026-10-16, identical commands)

Generated by `claude_patch_gmres6/analyse.py`.

The nonlinear solver, utils.py and the Firedrake installation have not changed since the 2026-10-16 star/single-column runs
(`data/paper_robustness/data_nonlinear_<test>_deg2_ew_{star,vanka}`), and the 440 commands are byte-identical (2 MPI
ranks). (Ten 2026-10-16 single-column dx runs at Nx = 80 started 17 s before the solver file's last edit that morning;
their logs and 9 of 10 counts match, so no effect is visible.) Differences between the two runs of a point are therefore
run-to-run noise. Both patches show rounding-level differences in every run with Nx >= 20 (none at Nx = 10), even
between the warm-up and the monitored pass inside one run. The noise check in ../60_mpi/README.md identifies the source:
repeats with a serial coarse LU (np = 1) or with SuperLU_DIST on the coarse level (np = 2) are bitwise identical, repeats
with the default parallel MUMPS coarse LU are not. The star's counts are insensitive to this noise; the single column's
slowly converging solves turn it into different counts. Paths relative to `yr2/`.

## Figures and tables
- `new_vs_old_counts.png`: new against old total FGMRES count per point, per patch, coloured by sweep.
- `repeat_counts.csv`: every point: new and old count/Newton steps or status, difference.

## Numbers (auto-generated)
| patch | points where both runs converged | identical counts | max abs difference | median abs difference | points where exactly one run failed | points where both failed with a different outcome (maxit vs diverged) |
|---|---|---|---|---|---|---|
''' + '\n'.join(lines) + '\n'
    open(os.path.join(out, 'README.md'), 'w').write(txt)
    print('wrote', out)


# --------------------------------------------------------------------------------------------- 5. MPI robustness
def mpi_folder():
    out = os.path.join(FIG, '60_mpi')
    os.makedirs(out, exist_ok=True)
    pts = list(csv.DictReader(open(os.path.join(MPI, 'points.csv'))))
    cache, res = {}, []
    for p in pts:
        cmd = p['command']
        g = lambda k, cast=float: cast(re.search(rf'--{k} (\S+)', cmd).group(1))   # noqa: E731
        test = re.search(r'--(dt|ar|dx|dz)_test', cmd).group(1)
        dt, nx, nz, h, shift = g('dt'), g('nx', int), g('nz', int), g('height'), g('shift')
        val = {'dt': shift * dt, 'ar': h, 'dx': nx, 'dz': nz}[test]
        d = os.path.join(MPI, p['folder'])
        if d not in cache:
            cache[d] = job_table(d)
        r = record(d, p['code'], test, val, dt, cache[d])
        res.append(dict(p, test=test, dt=dt, rec=r, np_=int(p['np'])))
    flat = [dict(code=x['code'], point=x['point'], description=x['description'], patch=x['patch'], np=x['np_'],
                 repeat=x['repeat'], count=fmt(x['rec']['count']), newton=fmt(x['rec']['newton']), status=x['rec']['status'],
                 secs=fmt(x['rec']['secs']), folder=x['folder']) for x in res]
    write_csv(os.path.join(out, 'mpi_counts.csv'), flat, list(flat[0].keys()))
    nps = sorted({x['np_'] for x in res})
    md = []
    for code in CODES:
        cp = [x for x in res if x['code'] == code]
        pids = sorted({x['point'] for x in cp}, key=lambda s: int(s[1:]))
        fig, axs = plt.subplots(2, 5, figsize=(16, 6.2), sharex=True)
        rows_md = []
        for ax, pid in zip(axs.flat, pids):
            desc = [x['description'] for x in cp if x['point'] == pid][0]
            for p in PATCHES:
                med, lo, hi, xs = [], [], [], []
                cells = []
                for n in nps:
                    v = [x['rec']['count'] for x in cp if x['point'] == pid and x['patch'] == p and x['np_'] == n and ok(x['rec'])]
                    st = [x['rec']['status'] for x in cp if x['point'] == pid and x['patch'] == p and x['np_'] == n]
                    if not st:
                        cells.append('-')
                        continue
                    if v:
                        xs.append(n); med.append(np.median(v)); lo.append(min(v)); hi.append(max(v))
                        cells.append(f'{int(np.median(v))}' + (f' [{min(v)}-{max(v)}]' if len(v) > 1 else '') +
                                     ('' if len(v) == len(st) else ' (' + ','.join(s for s in st if s != 'ok') + ')'))
                    else:
                        cells.append(','.join(st))
                if xs:
                    st_ = ps.STAR if p == 'star' else ps.COLUMN
                    ax.errorbar(xs, med, yerr=[np.array(med) - lo, np.array(hi) - med], capsize=3, label=PNAME[p],
                                color='C0' if p == 'star' else 'C3', **st_)
                # runs that ended without a count (aborted, time-out): a cross at the top of the panel
                for n in nps:
                    bad = [x['rec']['status'] for x in cp if x['point'] == pid and x['patch'] == p and x['np_'] == n
                           and x['rec']['status'] not in ('ok', 'maxit')]
                    if bad:
                        yb = 0.62 if p == 'star' else 0.50
                        ax.plot([n], [yb], 'x', ms=8, mew=2, color='C0' if p == 'star' else 'C3',
                                transform=ax.get_xaxis_transform(), clip_on=False)
                        ax.text(n, yb, f'{bad[0]}  ', ha='right', va='center', fontsize=7,
                                color='C0' if p == 'star' else 'C3', transform=ax.get_xaxis_transform())
                # reruns of aborted points with a working coarse solver (make_np16_check_jobs.py): open diamonds
                for v, (c, nn) in np16_check_counts(code, pid, p).items():
                    ax.plot([nn * (1.06 if 'superlu' in v else 0.94)], [c], 'D', ms=6, mfc='none', mew=1.3,
                            color='C0' if p == 'star' else 'C3',
                            label='np = 16 rerun (SuperLU_DIST or\nMUMPS ICNTL(14) = 200)' if (p == 'star' and 'superlu' in v) else None)
                rows_md.append([pid, desc, PNAME[p]] + cells)
            ax.set_xscale('log', base=2)
            ax.set_xticks(nps, [str(n) for n in nps])
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            ps.finish(ax, 'MPI ranks', 'iterations' if pid in ('P1', 'P6') else None,
                      legend=(pid == 'P1' or bool(np16_check_counts(code, pid, 'star'))), small=True, title=f'{pid}: {desc}')
            if not ax.lines:   # nothing converged at this point (e.g. linear P2): say what happened instead of an empty panel
                sts = sorted({x['rec']['status'] for x in cp if x['point'] == pid})
                ax.text(0.5, 0.5, 'every run of both patches, at every np:\n' + ', '.join(f'>= {MAXIT} iterations (maxit)' if t == 'maxit' else t for t in sts),
                        ha='center', va='center', fontsize=8, transform=ax.transAxes)
                ax.set_yticks([])
            else:
                lo_, hi_ = ax.get_ylim()
                if hi_ - lo_ < 3:   # flat lines: widen so that the integer ticks show
                    mid = 0.5 * (lo_ + hi_)
                    ax.set_ylim(mid - 2, mid + 2)
        for ax in list(axs.flat)[len(pids):]:
            ax.axis('off')
        fig.suptitle(f'{code} code: total outer iterations against MPI ranks (median, bars = min-max of repeats at np = 2 and 8)', fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'{code}_iterations_vs_np.png'), dpi=160)
        plt.close(fig)
        md.append(f'### {code} code\n\n' + md_table(['point', 'description', 'patch'] + [f'np={n}' for n in nps], rows_md))
        # error histories at P1 and P9 for every np (all repeats)
        fig, axs = plt.subplots(2, 2, figsize=(11, 6.5))
        for i, pid in enumerate(('P1', 'P9')):
            for j, p in enumerate(PATCHES):
                ax = axs[i, j]
                cols = ps.ramp(len(nps))
                for k, n in enumerate(nps):
                    for x in cp:
                        if x['point'] == pid and x['patch'] == p and x['np_'] == n and x['rec']['error'] is not None:
                            e = x['rec']['error']
                            ax.semilogy(np.arange(len(e)), np.maximum(e, 1e-16), color=cols[k], lw=1.0,
                                        label=f'np={n}' if x['repeat'] == '1' else None)
                desc = [x['description'] for x in cp if x['point'] == pid][0]
                ps.finish(ax, 'iteration', 'relative error', legend=(i == 0 and j == 0), small=True, title=f'{pid} {desc}, {PNAME[p]}')
        fig.tight_layout()
        fig.savefig(os.path.join(out, f'{code}_histories_P1_P9.png'), dpi=160)
        plt.close(fig)
    # pass criterion
    verdict = []
    for code in CODES:
        for pid in sorted({x['point'] for x in res if x['code'] == code}, key=lambda s: int(s[1:])):
            for p in PATCHES:
                xs = [x for x in res if x['code'] == code and x['point'] == pid and x['patch'] == p]
                st = {x['rec']['status'] for x in xs}
                if 'missing' in st:
                    verdict.append((code, pid, p, 'incomplete', '', '', ''))
                    continue
                if len(st) == 1 and st != {'ok'}:
                    verdict.append((code, pid, p, 'same status at every np', ','.join(st), '', ''))
                    continue
                # runs without a count (aborted etc.) are reported, the rest is judged on its own
                failed = sorted({(x['np_'], x['rec']['status']) for x in xs if not ok(x['rec'])})
                xs = [x for x in xs if ok(x['rec'])]
                if not xs:
                    verdict.append((code, pid, p, 'CHECK: no run converged: ' + ', '.join(f'np = {n}: {s}' for n, s in failed), '', '', ''))
                    continue
                # within-np spread: repeats at np = 2 and at np = 8; across-np spread: all runs of the point
                within = 0
                for n in (2, 8):
                    c = [x['rec']['count'] for x in xs if x['np_'] == n]
                    if len(c) > 1:
                        within = max(within, max(c) - min(c))
                allc = [x['rec']['count'] for x in xs]
                across = max(allc) - min(allc)
                v = 'no np effect' if across <= within + 2 else 'CHECK'
                if failed:
                    v = ', '.join(f'np = {n}: {s}' for n, s in failed) + f'; other runs: {v}'
                verdict.append((code, pid, p, v, f'{min(allc)}-{max(allc)}', f'{across}', f'{within}'))
    txt = f'''# MPI robustness: iteration counts against the number of MPI ranks

Generated by `claude_patch_gmres6/analyse.py` from `data/patch_gmres6/mpi/` (jobs made by `claude_patch_gmres6/make_mpi_jobs.py`,
one folder per code_patch / np<NN>_r<repeat>; the exact command of every run is in `data/patch_gmres6/mpi/points.csv` and in
each `jobs/job_*/cmd.sh`).

Same solver configurations as the main sweeps (see ../11_linear_dt/README.md and ../21_nonlinear_dt/README.md), on
np = {nps} MPI ranks (np <= Nx, the number of coarse cells: Nx = 10 runs only up to np = 8). Points P1-P9 (per code, see the
tables) are the nominal case, the extremes of each sweep and one known single-column stress point. np = 2 and np = 8 were
run three times each to measure the run-to-run noise (the parallel MUMPS coarse LU is not bitwise reproducible; np = 1 uses
PETSc's serial LU and is deterministic).

Why counts can differ with np at all: the patches (and hence the preconditioner) are the same for every np in exact
arithmetic (Firedrake builds one patch per owned base vertex/cell; verified at np = 2, 4), but the floating-point order of
operations changes with the partition (and the coarse LU is PETSc's serial LU at np = 1, MUMPS at np >= 2); with a Krylov
smoother inside FGMRES the preconditioner is not a fixed linear operator, and where an iteration converges slowly such
rounding differences can change the count (see ../40_repeatability_nonlinear/).

Criterion used below: for each (code, point, patch), every run must end with the same status, and the spread of the
counts over ALL runs (every np) must not exceed the largest spread between repeats at a fixed np (np = 2 or np = 8) by
more than 2 iterations, i.e. changing the number of ranks must not move the count more than repeating the run does.

## Figures and tables
- `<code>_iterations_vs_np.png`: one panel per point: median count against np, bars = min-max over the repeats (np = 2, 8).
- `<code>_histories_P1_P9.png`: relative-error histories of every run at the nominal point P1 and the stress point P9.
- `mpi_counts.csv`: every run (code, point, patch, np, repeat, count, Newton steps, status, run time, folder).

## Verdict (auto-generated)
{md_table(['code', 'point', 'patch', 'verdict', 'counts over all runs (or statuses)', 'spread over all np', 'spread between repeats at fixed np'], verdict)}

## Counts: median [min-max over repeats] per np
''' + '\n\n'.join(md) + '\n' + np16_check_md() + noise_check_md()
    open(os.path.join(out, 'README.md'), 'w').write(txt)
    print('wrote', out)


def np16_check_runs():
    '''Reruns of the aborted np = 16 runs (make_np16_check_jobs.py): {(patch, variant): (record, run.log text)}.'''
    out = {}
    for d in sorted(glob.glob(os.path.join(MPI, 'np16_check', '*_*'))):
        patch, variant = os.path.basename(d).split('_', 1)
        logs = glob.glob(os.path.join(d, 'jobs', '*', 'run.log'))
        out[(patch, variant)] = (record(d, 'nonlinear', 'dx', 80, 50.0), open(logs[0]).read() if logs else '')
    return out


def np16_check_counts(code, pid, patch):
    '''{variant: (count, np)} of the reruns that converged, for the panel of (code, pid, patch).'''
    if code != 'nonlinear' or pid != 'P6':
        return {}
    return {v: (r['count'], 16) for (p, v), (r, _) in np16_check_runs().items() if p == patch and ok(r)}


def np16_check_md():
    runs = np16_check_runs()
    if not runs:
        return ''
    label = {'mumps': 'MUMPS (default), converged-reason views on',
             'superlu': 'SuperLU_DIST on the coarse level',
             'mumps_icntl14': 'MUMPS with ICNTL(14) = 200 (200 % extra working space; default 20 %)'}
    rows = []
    for (patch, variant), (r, log) in sorted(runs.items(), key=lambda kv: (kv[0][0], list(label).index(kv[0][1]))):
        why = re.findall(r'Linear (\S+) solve did not converge due to (\S+)[^\n]*\n\s*PC failed due to (\S+)', log)
        first = f'{why[0][0]}: {why[0][1]}, PC failed due to {why[0][2]}' if why else ''
        rows.append([PNAME[patch], label[variant], cell(r) if ok(r) else 'aborted', first])
    return ('\n## The two aborted runs (nonlinear code, P6 = Nx 80, np = 16)\n\n'
            'Both patches stopped in the first Newton step with a LAPACK error inside the ksp_converged_rate diagnostic '
            '(the outer KSP had returned after 0 iterations). Reruns (`make_np16_check_jobs.py`, data in '
            '`data/patch_gmres6/mpi/np16_check/`) with converged-reason views show why: the MUMPS LU factorisation of '
            'the multigrid coarse level fails (FACTOR_OUTMEMORY: MUMPS ran out of the working space it had reserved; '
            'the machine had more than 30 GB of free memory), the coarse solve returns NaN/Inf and the whole '
            'preconditioner fails. With SuperLU_DIST on the coarse level, or with MUMPS and 200 % extra working space, '
            'both patches converge at np = 16 with the same counts (and Newton steps) as at np = 1-8. The linear code '
            'at the same point and every other point at np = 16 ran with the default settings. So this is a '
            'coarse-solver setting, not an effect of the patch or of the column smoother.\n\n'
            + md_table(['patch', 'coarse solver', 'count (total FGMRES / Newton steps)', 'first failing solve'], rows) + '\n')


def noise_check_md():
    '''Counts and bitwise equality of the output files across repeats of two nonlinear single-column points (P3, P9):
    serial with PETSc's LU on the coarse level, 2 ranks with MUMPS (the default), 2 ranks with SuperLU_DIST.'''
    nc = os.path.join(MPI, 'noise_check')
    if not os.path.isdir(nc):
        return ''
    groups = {
        'np = 1, PETSc LU (coarse)': [os.path.join(MPI, 'nonlinear_vanka', 'np01_r1')] + sorted(glob.glob(os.path.join(nc, 'np01_petscLU_r*'))),
        'np = 2, MUMPS (default)': sorted(glob.glob(os.path.join(MPI, 'nonlinear_vanka', 'np02_r*'))),
        'np = 2, SuperLU_DIST (coarse)': sorted(glob.glob(os.path.join(nc, 'np02_superlu_r*'))),
    }
    pts = {'P3': ('dt', 1e-5, 50000.0, 'shift'), 'P9': ('dz', 80, 5000.0, 'nz')}
    rows = []
    for pid, (test, val, dt, key) in pts.items():
        for gname, ds in groups.items():
            counts, blobs = [], []
            for d in ds:
                r = record(d, 'nonlinear', test, val, dt)
                counts.append(cell(r))
                f = find(d, 'residual', dt, key, target(test, val, dt))
                blobs.append(open(f, 'rb').read() if f else None)
            done = [b for b in blobs if b is not None]
            ident = ('yes' if len(set(done)) == 1 else 'no') if len(done) > 1 else 'n/a'
            rows.append([pid, gname, ', '.join(counts), f'{len(done)} runs: {ident}'])
    # did the SuperLU option take effect? (PETSc lists unused options at the end of run.log)
    unused = []
    for log in glob.glob(os.path.join(nc, 'np02_superlu_r*', 'jobs', '*', 'run.log')):
        t = open(log).read()
        if 'superlu' in t.split('options you set that were not used')[-1] and 'not used' in t:
            unused.append(log)
    note = ('SuperLU_DIST option reported as UNUSED by PETSc in: ' + ', '.join(unused)) if unused else \
        'PETSc reported no unused option in the SuperLU_DIST runs (the coarse-solver override took effect).'
    return ('\n## Noise check: what makes the single-column counts vary?\n\n'
            'Two nonlinear single-column points with K < 1 (P3: dt = 5e4 s, C1 = 1e-5; P9: Nz = 80, dt = 5000 s), repeated '
            '(jobs from `make_noise_jobs.py`, data in `data/patch_gmres6/mpi/noise_check/` and the MPI-study folders). '
            '"bitwise identical" compares the residual_*.out files (full-precision residual history) of the repeats. ' + note + '\n\n'
            + md_table(['point', 'runs', 'counts (total FGMRES / Newton steps)', 'bitwise identical'], rows) + '\n')


# --------------------------------------------------------------------------------------------- 6. overview
def overview_folder(allrows):
    out = os.path.join(FIG, '00_overview')
    os.makedirs(out, exist_ok=True)
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.8))
    lines = []
    for ax, code in zip(axs, CODES):
        data = []
        for test in TESTS:
            rows = allrows[(code, test)]
            v = [ratio(r['vanka'], r['star']) for r in rows]
            data.append([x for x in v if np.isfinite(x)] or [np.nan])
            fs = sum(1 for r in rows if r['star']['status'] not in ('ok', 'missing'))
            fv = sum(1 for r in rows if r['vanka']['status'] not in ('ok', 'missing'))
            ms = sum(1 for r in rows if r['star']['status'] == 'missing') + sum(1 for r in rows if r['vanka']['status'] == 'missing')
            lines.append(f'| {code} | {test} | {stats_line(v)} | {fs} | {fv} | {ms} |')
        ax.boxplot(data, tick_labels=list(TESTS))
        ax.axhline(1.0, color=ps.GREY, lw=0.8, ls=':')
        ps.finish(ax, 'sweep', 'single column / star', legend=False, title=f'{code} code, GMRES(6): iteration ratio per point')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'ratio_by_sweep.png'), dpi=160)
    plt.close(fig)
    # ratio against K and rho1 (all points, both codes)
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
    npts = 0
    tcol = {'dt': 'C0', 'ar': 'C1', 'dx': 'C2', 'dz': 'C3'}
    for code, mk in (('linear', 'o'), ('nonlinear', 's')):
        for t in TESTS:
            rows = [r for r in allrows[(code, t)] if np.isfinite(ratio(r['vanka'], r['star']))]
            npts += len(rows)
            if not rows:
                continue
            K = [r['K'] for r in rows]; rho = [r['rho1'] for r in rows]; q = [ratio(r['vanka'], r['star']) for r in rows]
            kw = dict(s=16, marker=mk, alpha=0.75, label=f'{code} {t}')
            kw.update(dict(facecolors='none', edgecolors=tcol[t]) if code == 'linear' else dict(color=tcol[t]))
            axs[0].scatter(K, q, **kw)
            axs[1].scatter(rho, q, **kw)
    for ax, xl in zip(axs, (r'$K = \delta / \delta_{thr}$', r'$\rho_1 = m L_z^2 / (\pi^2 \Delta x^2)$')):
        if npts:
            ax.set_xscale('log')
        ax.axhline(1.0, color=ps.GREY, lw=0.8, ls=':')
        ps.finish(ax, xl, 'single column / star', small=True, legend=(ax is axs[0]), ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'ratio_vs_K_rho1.png'), dpi=160)
    plt.close(fig)
    txt = '''# Overview: star vs single-column patch with GMRES(6) smoothing (2026-11-16)

Generated by `claude_patch_gmres6/analyse.py`. One row per (code, sweep). Detailed folders:
`11-14_linear_<test>`, `21-24_nonlinear_<test>` (comparison per sweep), `30_smoother_linear` (GMRES(6) vs the earlier
Chebyshev(6) linear runs), `40_repeatability_nonlinear` (run-to-run noise), `50_paper_style` (the paper's ten figures for
each patch), `60_mpi` (MPI ranks). The study write-up is `../README.md`.

## Figures
- `ratio_by_sweep.png`: distribution of single column / star over all points of each sweep (both converged).
- `ratio_vs_K_rho1.png`: the same ratios against K = delta/delta_thr and rho1 = m Lz^2/(pi^2 dx^2); colour = sweep, open
  circles = linear code, filled squares = nonlinear code.

Notes. (1) The nominal configuration (Nx = Nz = 20, Lz = 4 km, the sweep's default C1) is a member of all four sweeps, so
pooled statistics over sweeps count it up to four times. (2) K is computed with the same formula for both codes; the
earlier study (`Boussinesq_implicit/docs/claude_studies/claude_small_patch_robustness.md`, end of Section 4.4) derives,
without measuring it, that because of the Irksome stage scaling the nonlinear script's delta is worth 2 delta of the linear
code (K_eff = 2K), so like-for-like comparisons between the codes may be off by a factor of 2 in K.

## Numbers (auto-generated)
| code | sweep | single column / star | star failures | single-column failures | not (yet) run |
|---|---|---|---|---|---|
''' + '\n'.join(lines) + '\n'
    open(os.path.join(out, 'README.md'), 'w').write(txt)
    print('wrote', out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', nargs='+', default=['main', 'paper', 'smoother', 'repeat', 'mpi', 'overview'])
    a = ap.parse_args()
    ps.apply()
    allrows = {}
    if {'main', 'smoother', 'repeat', 'overview'} & set(a.only):
        for code in CODES:
            for test in TESTS:
                allrows[(code, test)] = main_folder(code, test) if 'main' in a.only else collect(code, test)
    if 'paper' in a.only:
        paper_folder()
    if 'smoother' in a.only:
        smoother_folder(allrows)
    if 'repeat' in a.only:
        repeat_folder(allrows)
    if 'mpi' in a.only and os.path.exists(os.path.join(MPI, 'points.csv')):
        mpi_folder()
    if 'overview' in a.only:
        overview_folder(allrows)
