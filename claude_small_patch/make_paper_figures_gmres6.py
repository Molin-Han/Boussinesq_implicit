'''
The paper's iteration figures and counts from the 2026-11-16 re-runs (GMRES(6) level smoother everywhere, pressure
back-substitution in the nonlinear solver, nonlinear shift delta = C_1/Delta t in the paper's convention), in one
counting convention:

  * linear runs:    n = number of outer FGMRES iterations to reach the stopping test (the error file has n + 1 entries,
                    iterations 0..n); histories show every iterate, the last one included;
  * nonlinear runs: T = total number of Krylov iterations over the Newton steps of the monitored time step
                    (snes_ksp_cum[-1]; Newton steps are not added); histories show every Newton iterate;
  * the error is measured against the converged solution of an identical warm-up solve, so the last point of every
    history is at rounding level; it is drawn at max(error, 1e-16);
  * count plots omit runs in which a Krylov solve reached the cap of 100 iterations or that did not finish.

Data: <data>/data_linear_<test>_deg2[_vanka], <data>/data_nonlinear_<test>_deg2_ew[_vanka] (made by the run_one_*.sh
scripts there), default /Users/molin_han/Code_Repo/data/paper_robustness_gmres6.  The nonlinear code takes its shift in
the stage form of the implicit midpoint rule, where it acts as twice the shift of A_delta; the runs therefore used
--shift C_1/(2 Delta t), and the files of the nonlinear shift sweep carry that value in their names.

Writes (same names and sizes as the figures the paper includes, so the tex does not change):
  <paper>/Image_restyled/half/error_{dx,dz,AR}_Robust_dt50.0.png, half/error_shift.png        (fig:linear)
  <paper>/Image_restyled/snes_error_{dx,dz}_dt50.0.png, snes_error_AR_dt100.0.png,
  <paper>/Image_restyled/error_shift_nonlinear.png                                              (fig:nonlinear)
  <paper>/Image_restyled/single_column_{linear,nonlinear}.png                                   (supplement)
  <data>/counts_gmres6.csv                                                                      (every run)

    python make_paper_figures_gmres6.py [--data DIR] [--paper DIR]
'''
import argparse
import csv
import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_style as ps  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

LX = 3.0e5
MAXIT = 100
FLOOR = 1e-16
LIN = dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000.], C1=1e-3, C1_list=[1., 0.1, 0.01, 1e-3, 1e-4, 1e-5],
           heights=[12000, 8000, 4000, 2000, 1000, 500, 250], nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80])
NL = dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000., 50000.], C1=1e-5,
          C1_list=[1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6], heights=[12000, 8000, 4000, 2000, 1000],
          nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80])
SWEEP = {'linear': LIN, 'nonlinear': NL}


def folder(data, mode, test, patch):
    return f'{data}/data_{mode}_{test}_deg2' + ('_ew' if mode == 'nonlinear' else '') + ('_vanka' if patch == 'vanka' else '')


def find(d, prefix, dt, key, tgt):
    '''file <prefix>_dt<dt>_<key><value>.out whose value is closest (in log) to tgt, within 20 %.'''
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


def load(d, prefix, dt, key, tgt):
    f = find(d, prefix, dt, key, tgt)
    return None if f is None else np.atleast_1d(np.loadtxt(f))


def members(mode, test, short=True):
    '''(legend label, file key, sweep value, file-name target) per member.'''
    S = SWEEP[mode]
    if test == 'dx':
        return [(rf'$N_x = {n}$', 'nx', n, n) for n in S['nxs']]
    if test == 'dz':
        return [(rf'$N_z = {n}$', 'nz', n, n) for n in S['nzs']]
    if test == 'ar':
        return [(rf'$L_z = {h/1000:g}$ km', 'ar', h / LX, h / LX) for h in S['heights']]
    return [(r'$C_1 = 1$' if c == 1 else rf'$C_1 = {ps.sci(c)[1:-1]}$', 'shift', c, None) for c in S['C1_list']]


def shift_target(mode, c1, dt):
    # the nonlinear runs pass --shift C_1/(2 dt) (stage form), the linear ones --shift C_1/dt
    return c1 / dt if mode == 'linear' else c1 / (2.0 * dt)


def job_index(d):
    '''rc and converged reasons of every job of a data folder, keyed by the solver arguments.'''
    idx = {}
    for cmd in glob.glob(f'{d}/jobs/*/cmd.txt'):
        jd = os.path.dirname(cmd)
        args = open(cmd).read().split()
        rc = open(f'{jd}/rc').read().strip() if os.path.exists(f'{jd}/rc') else None
        a = {args[i][2:]: args[i + 1] for i in range(1, len(args) - 1) if args[i].startswith('--') and not args[i + 1].startswith('--')}
        idx[jd] = dict(args=a, rc=rc, log=f'{jd}/run.log')
    return idx


def job_for(idx, test, key, val, dt, mode, c1=None):
    for jd, j in idx.items():
        a = j['args']
        if float(a.get('dt', 'nan')) != dt:
            continue
        if test == 'dt':
            if abs(np.log(float(a['shift'])) - np.log(shift_target(mode, c1, dt))) < 1e-6:
                return j
        elif test == 'ar':
            if abs(float(a['height']) / float(a['length']) - val) < 1e-9:
                return j
        elif test == 'dx':
            if int(a['nx']) == val:
                return j
        elif test == 'dz':
            if int(a['nz']) == val:
                return j
    return None


def linear_reason(log):
    '''KSP converged reason of the monitored (last) linear solve, as printed by LB_time_slice.py.'''
    if not os.path.exists(log):
        return None
    r = None
    for line in open(log, errors='replace'):
        m = re.match(r'^\*+\s+(-?\d+)\s*$', line.strip())
        if m:
            r = int(m.group(1))
    return r


def record(data, mode, test, patch, key, val, tgt, dt, idx):
    d = folder(data, mode, test, patch)
    rec = dict(mode=mode, test=test, patch=patch, member=val, dt=dt, count=np.nan, newton=np.nan, status='missing')
    c1 = val if test == 'dt' else SWEEP[mode]['C1']
    t = shift_target(mode, c1, dt) if test == 'dt' else tgt
    job = job_for(idx, test, key, val, dt, mode, c1)
    if mode == 'linear':
        e = load(d, 'error', dt, key, t)
        if e is not None:
            n = len(e) - 1
            reason = linear_reason(job['log']) if job else None
            if reason is not None and reason < 0:
                rec['status'] = 'div'
            elif n >= MAXIT or reason == 4:
                rec['status'] = 'maxit'
            else:
                rec['status'], rec['count'] = 'ok', n
    else:
        cu = load(d, 'snes_ksp_cum', dt, key, t)
        if cu is not None and cu.size > 0:
            per = np.diff(cu) if cu.size > 1 else cu
            rec['newton'] = cu.size - 1
            if np.any(per >= MAXIT):
                rec['status'] = 'maxit'
            else:
                rec['status'], rec['count'] = 'ok', int(cu[-1])
    if rec['status'] == 'missing' and job is not None and job['rc'] not in (None, '0'):
        rec['status'] = 'div'
    return rec


def all_records(data):
    recs = []
    for mode in ('linear', 'nonlinear'):
        S = SWEEP[mode]
        for test in ('dt', 'ar', 'dx', 'dz'):
            for patch in ('star', 'vanka'):
                idx = job_index(folder(data, mode, test, patch))
                for (_, key, val, tgt) in members(mode, test):
                    for dt in S['dts']:
                        recs.append(record(data, mode, test, patch, key, val, tgt, dt, idx))
    return recs


def lookup(recs, mode, test, patch, val, dt):
    for r in recs:
        if (r['mode'], r['test'], r['patch'], r['dt']) == (mode, test, patch, dt) and np.isclose(r['member'], val, rtol=1e-9):
            return r
    raise KeyError((mode, test, patch, val, dt))


# ---------------------------------------------------------------------------------------------------------------
def linear_history(data, test, out):
    fig, ax = plt.subplots(figsize=ps.HALF)
    mem = members('linear', test)
    d = folder(data, 'linear', test, 'star')
    for i, ((lab, key, val, tgt), col) in enumerate(zip(mem, ps.ramp(len(mem)))):
        e = load(d, 'error', 50.0, key, tgt)
        if e is None:
            continue
        if test == 'ar':   # seven members: label the thickest, a middle and the thinnest domain (colour dark to light)
            lab = f'{val * LX / 1000:g}' if i in (0, len(mem) // 2, len(mem) - 1) else None
        ax.semilogy(np.arange(len(e)), np.clip(e, FLOOR, None), color=col, label=lab, **ps.STAR)
    ps.finish(ax, xlabel='FGMRES iteration', ylabel='relative error', legend=(test != 'ar'))
    if test == 'ar':   # in the empty space between the drops of the 4 km and 8 km curves
        ax.legend(title=r'$L_z$ (km)', loc='lower left', bbox_to_anchor=(0.43, 0.10), borderaxespad=0.0, ncol=1, fontsize=7,
                  title_fontsize=7, handlelength=1.0, handletextpad=0.4, borderpad=0.3, labelspacing=0.2)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.savefig(out); plt.close(fig)


def linear_shift_sweep(recs, out):
    fig, ax = plt.subplots(figsize=ps.HALF)
    mem = members('linear', 'dt')
    for (lab, _, c1, _), col in zip(mem, ps.ramp(len(mem))):
        xs = [dt for dt in LIN['dts'] if lookup(recs, 'linear', 'dt', 'star', c1, dt)['status'] == 'ok']
        ax.plot(xs, [lookup(recs, 'linear', 'dt', 'star', c1, dt)['count'] for dt in xs], color=col, label=lab, **ps.STAR)
    ax.set_xscale('log')
    ps.finish(ax, xlabel=r'$\Delta t$ (s)', ylabel='FGMRES iterations', small=True)
    ax.legend(loc='upper right', ncol=2, columnspacing=0.6, handlelength=1.2, borderaxespad=0.2)
    top = max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata()))
    ax.set_ylim(0, 1.65 * top)
    fig.savefig(out); plt.close(fig)


def nonlinear_history(data, test, dt, out):
    fig, ax = plt.subplots(figsize=ps.HALF)
    mem = members('nonlinear', test)
    d = folder(data, 'nonlinear', test, 'star')
    for (lab, key, val, tgt), col in zip(mem, ps.ramp(len(mem))):
        se = load(d, 'snes_error', dt, key, tgt)
        cu = load(d, 'snes_ksp_cum', dt, key, tgt)
        if se is None or cu is None or se.size == 0:
            continue
        ax.semilogy(cu, np.clip(se, FLOOR, None), color=col, label=lab, **ps.STAR)
    ps.finish(ax, xlabel='cumulative FGMRES iterations', ylabel=r'$\|U_k - U^*\| / \|U^*\|$', legend=False)
    ax.legend(loc='best', **LEG8)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.savefig(out); plt.close(fig)


def nonlinear_shift_sweep(recs, out):
    fig, ax = plt.subplots(figsize=ps.HALF)
    mem = members('nonlinear', 'dt')
    for (lab, _, c1, _), col in zip(mem, ps.ramp(len(mem))):
        xs = [dt for dt in NL['dts'] if lookup(recs, 'nonlinear', 'dt', 'star', c1, dt)['status'] == 'ok']
        ax.plot(xs, [lookup(recs, 'nonlinear', 'dt', 'star', c1, dt)['count'] for dt in xs], color=col, label=lab, **ps.STAR)
    ax.set_xscale('log')
    ps.finish(ax, xlabel=r'$\Delta t$ (s)', ylabel='total FGMRES iterations', legend=False)
    ax.legend(loc='upper right', **LEG8)
    top = max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata()))
    ax.set_ylim(0, 1.45 * top)   # room for the seven-entry legend above the curves right of the peak
    fig.savefig(out); plt.close(fig)


LEG = dict(fontsize=7.5, handlelength=1.6, borderpad=0.3, labelspacing=0.25, columnspacing=0.8)
LEG8 = dict(ps.LEGEND_SMALL, fontsize=8)   # the half-width nonlinear panels: 8 pt, so that sub- and superscripts print legibly


def single_column(recs, mode, out):
    S = SWEEP[mode]
    ylab = 'FGMRES iterations' if mode == 'linear' else 'total FGMRES iterations'
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.6))

    def sweep_panel(ax, test, title):
        mem = members(mode, test)
        cols = ps.ramp(len(mem))
        show = {0, len(mem) // 2, len(mem) - 1}
        for i, (lab, _, val, _) in enumerate(mem):
            for patch, style in (('star', ps.STAR), ('vanka', ps.COLUMN)):
                xs = [dt for dt in S['dts'] if lookup(recs, mode, test, patch, val, dt)['status'] == 'ok']
                ax.plot(xs, [lookup(recs, mode, test, patch, val, dt)['count'] for dt in xs], color=cols[i],
                        label=lab if (patch == 'star' and i in show) else None, **style)
        ax.set_xscale('log'); ax.set_title(title)
        ax.set_xlabel(r'$\Delta t$ (s)')

    ax = axes[0]
    sweep_panel(ax, 'ar', f"(a) aspect ratio, $C_1 = {ps.sci(S['C1'])[1:-1]}$")
    ps.patch_legend(ax, star_label='star', column_label='single column')
    ax.set_ylim(0, 1.75 * max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(ncol=1, loc='upper right', **LEG)
    ax.set_ylabel(ylab)
    ax = axes[1]
    for test, mk, colr, name in (('dx', 'o', 'C0', '$N_x$'), ('dz', 's', 'C1', '$N_z$')):
        ns = S['nxs'] if test == 'dx' else S['nzs']
        for patch, ls, mfc, pname in (('star', '-', colr, 'star'), ('vanka', '--', 'white', 'single column')):
            ys = [lookup(recs, mode, test, patch, n, 50.0)['count'] for n in ns]
            ax.plot(ns, ys, ls=ls, marker=mk, color=colr, mfc=mfc, label=f'{pname}, {name}')
    ax.set_title(r"(b) mesh resolution, $\Delta t = 50$ s")
    ax.set_xscale('log'); ax.set_xticks([10, 20, 40, 60, 80]); ax.set_xticklabels(['10', '20', '40', '60', '80'])
    ax.xaxis.set_minor_locator(plt.NullLocator()); ax.set_xlabel(r'$N_x$ (coarsest mesh) or $N_z$')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(0, 1.6 * max(np.nanmax(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(loc='lower right', ncol=1, **LEG)
    ax = axes[2]
    sweep_panel(ax, 'dt', r"(c) shift $\delta = C_1/\Delta t$")
    ax.set_ylim(0, 1.65 * max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(ncol=1, loc='upper right', **LEG)
    for a_ in axes:
        a_.grid(True, alpha=0.3)
    fig.tight_layout(w_pad=0.8)
    fig.savefig(out); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default='/Users/molin_han/Code_Repo/data/paper_robustness_gmres6')
    ap.add_argument('--paper', default=os.path.abspath(os.path.join(HERE, '..', '..', '..', 'Tex_Repo', 'Patch-Preconditioner-for-Boussinesq')))
    ap.add_argument('--out', default=None, help='write all figures into this one folder instead of the paper tree')
    a = ap.parse_args()
    img = os.path.join(a.paper, 'Image_restyled')
    path = (lambda sub, name: os.path.join(a.out, name)) if a.out else (lambda sub, name: os.path.join(img, sub, name))
    if a.out:
        os.makedirs(a.out, exist_ok=True)
    ps.apply()
    recs = all_records(a.data)
    with open(os.path.join(a.data, 'counts_gmres6.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['mode', 'test', 'patch', 'member', 'dt', 'count', 'newton_steps', 'status'])
        for r in recs:
            w.writerow([r['mode'], r['test'], r['patch'], r['member'], r['dt'], r['count'], r['newton'], r['status']])
    for test, name in (('dx', 'error_dx_Robust_dt50.0.png'), ('dz', 'error_dz_Robust_dt50.0.png'), ('ar', 'error_AR_Robust_dt50.0.png')):
        linear_history(a.data, test, path('half', name))
    linear_shift_sweep(recs, path('half', 'error_shift.png'))
    nonlinear_history(a.data, 'dx', 50.0, path('', 'snes_error_dx_dt50.0.png'))
    nonlinear_history(a.data, 'dz', 50.0, path('', 'snes_error_dz_dt50.0.png'))
    nonlinear_history(a.data, 'ar', 100.0, path('', 'snes_error_AR_dt100.0.png'))
    nonlinear_shift_sweep(recs, path('', 'error_shift_nonlinear.png'))
    single_column(recs, 'linear', path('', 'single_column_linear.png'))
    single_column(recs, 'nonlinear', path('', 'single_column_nonlinear.png'))
    print('figures written', a.out or img, '| counts', os.path.join(a.data, 'counts_gmres6.csv'))


if __name__ == '__main__':
    main()
