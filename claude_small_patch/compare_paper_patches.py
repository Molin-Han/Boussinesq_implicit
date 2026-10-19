'''
Item 1: the paper's robustness sweeps (delta = C_1/dt) with the star and the single-column (vanka) patch.

Reads
    ../data_<mode>_<test>_deg2[_ew]_{star,vanka}     the new sweeps (run_all_sweeps.sh / run_claude_paper.sh -g Bash/bash_paper_patch.py)
    ../data_<mode>_<test>_deg2[_ew]                  the paper's own data (star patch, 2026-05/06)
for mode in {linear, nonlinear}, test in {dt, ar, dx, dz}, and writes into results/:
    summary_<mode>_<test>.md            iteration counts: paper | star (re-run) | single column, per dt and member
    counts_<mode>_<test>.png            counts against dt, star (solid) and single column (dashed)
    history_<mode>_<test>_dt<dt>.png    error histories at the paper's dt (the figures of the paper, both patches)
    overview_<mode>.png                 four panels in the style of the paper's fig. 14, for the paper's parametrisation
    all_counts.csv                      every run of every folder

The figures show only the new runs (star re-run and single column); the paper's own data enter the tables only.

Counting conventions: linear n = number of KSP iterations to the stopping test (= lines of the error file - 1;
the paper's plots show the last plotted index n - 1); nonlinear = total KSP iterations over all Newton steps
(snes_ksp_cum[-1]); 'maxit' = a KSP solve of the run hit --maxit 100; 'div' = Newton diverged (no output file,
job status ConvergenceError); 'missing' = no file and no failed job (not run yet).

    python compare_paper_patches.py [--modes linear nonlinear] [--tests dt ar dx dz]
'''
import os
import re
import glob
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_style as ps  # noqa: E402
ps.apply()
TITLES = False   # set by --titles: working titles on the figures (off for the paper)
YR2 = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.join(HERE, 'results')
os.makedirs(OUT, exist_ok=True)
N_FREQ, LX, LZ = 1.0e-2, 3.0e5, 4.0e3

# the paper's sweeps (Bash/bash.py, Bash/bash_nonlinear.py)
SWEEP = {
    'linear': dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000.], C1=1e-3,
                   C1_list=[1., 0.1, 0.01, 1e-3, 1e-4, 1e-5], heights=[12000, 8000, 4000, 2000, 1000, 500, 250],
                   nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80], paper_dt={'ar': [50.], 'dx': [50.], 'dz': [50.], 'dt': []}),
    'nonlinear': dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000., 50000.], C1=1e-5,
                      C1_list=[1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6], heights=[12000, 8000, 4000, 2000, 1000],
                      nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80], paper_dt={'ar': [100., 5000.], 'dx': [50.], 'dz': [50.], 'dt': [500.]}),
}
MAXIT = 100


def folder(mode, test, patch):
    base = f'{YR2}/data_{mode}_{test}_deg2' + ('_ew' if mode == 'nonlinear' else '')
    return base if patch == 'paper' else f'{base}_{patch}'


def members(mode, test):
    S = SWEEP[mode]
    if test == 'dt':
        return [(rf'$C_1 = {ps.sci(c)[1:-1]}$', 'shift', c) for c in S['C1_list']]
    if test == 'ar':
        return [(rf'$L_z = {h/1000:g}$ km', 'ar', h) for h in S['heights']]
    if test == 'dx':
        return [(rf'$N_x = {n}$', 'nx', n) for n in S['nxs']]
    return [(rf'$N_z = {n}$', 'nz', n) for n in S['nzs']]


def target(test, val, dt):
    if test == 'dt':
        return val / dt            # shift = C1/dt
    if test == 'ar':
        return val / LX            # ar = height/length
    return val


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


def failed_jobs(d):
    out = []
    for st in glob.glob(f'{d}/jobs/*/status'):
        try:
            rc, secs, reason = open(st).read().split()[:3]
        except ValueError:
            continue
        if rc == '0':
            continue
        cmd = open(os.path.join(os.path.dirname(st), 'cmd.sh')).read()
        g = lambda k, cast=float: cast(re.search(rf'--{k} (\S+)', cmd).group(1))
        out.append(dict(dt=g('dt'), shift=g('shift'), height=g('height'), nx=g('nx', int), nz=g('nz', int), reason=reason))
    return out


def record(mode, test, patch, key, val, dt, failed):
    d = folder(mode, test, patch)
    tgt = target(test, val, dt)
    rec = dict(mode=mode, test=test, patch=patch, member=val, dt=dt, count=np.nan, newton=np.nan, status='missing',
               error=None, snes_error=None, snes_cum=None)
    err = find(d, 'error', dt, key, tgt)
    if err is not None:
        rec['error'] = np.atleast_1d(np.loadtxt(err))
    if mode == 'linear':
        if err is not None:
            n = len(rec['error']) - 1
            rec['status'] = 'maxit' if n >= MAXIT else 'ok'
            rec['count'] = n if n < MAXIT else np.nan
    else:
        cumf = find(d, 'snes_ksp_cum', dt, key, tgt)
        if cumf is not None:
            cum = np.atleast_1d(np.loadtxt(cumf)).astype(int)
            rec['snes_cum'] = cum
            se = find(d, 'snes_error', dt, key, tgt)
            rec['snes_error'] = np.atleast_1d(np.loadtxt(se)) if se else None
            rec['newton'] = len(cum) - 1
            per = np.diff(cum) if len(cum) > 1 else cum
            rec['status'] = 'maxit' if np.any(per >= MAXIT) else 'ok'
            rec['count'] = int(cum[-1]) if rec['status'] == 'ok' else np.nan
    if rec['status'] == 'missing':
        for f in failed:
            if abs(f['dt'] - dt) > 1e-9 * max(1., dt):
                continue
            if key == 'shift' and abs(np.log(f['shift']) - np.log(tgt)) > 0.2:
                continue
            if key == 'ar' and abs(f['height'] - val) > 1e-6:
                continue
            if key == 'nx' and f['nx'] != val:
                continue
            if key == 'nz' and f['nz'] != val:
                continue
            rec['status'] = 'div' if ('Converg' in f['reason'] or 'DIVERGED' in f['reason']) else f['reason']
    return rec


def cell(r):
    if r['status'] == 'ok':
        return f"{int(r['count'])}" + (f"/{int(r['newton'])}" if np.isfinite(r['newton']) else '')
    return r['status']


def main(modes, tests):
    allrows = []
    for mode in modes:
        S = SWEEP[mode]
        panels = {}
        for test in tests:
            mem = members(mode, test)
            data = {}
            for patch in ['paper', 'star', 'vanka']:
                d = folder(mode, test, patch)
                if not os.path.isdir(d):
                    continue
                fj = failed_jobs(d)
                data[patch] = {(val, dt): record(mode, test, patch, key, val, dt, fj) for (_, key, val) in mem for dt in S['dts']}
                for r in data[patch].values():
                    allrows.append(r)
            if not data:
                continue
            # ---- markdown table
            lines = [f'# {mode} {test} test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column\n',
                     ('linear: KSP iterations to rtol 1e-8' if mode == 'linear' else 'nonlinear: total KSP iterations over all Newton steps / Newton steps; f = 0, Eisenstat-Walker on')
                     + f". Folders: {', '.join(os.path.basename(folder(mode, test, p)) for p in data)}.\n",
                     '| dt | ' + ' | '.join(lab for (lab, _, _) in mem) + ' |', '|---|' + '---|' * len(mem)]
            for dt in S['dts']:
                cells = []
                for (_, _, val) in mem:
                    cells.append(' / '.join(cell(data[p][(val, dt)]) for p in ['paper', 'star', 'vanka'] if p in data))
                lines.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')
            # ratios
            if 'star' in data and 'vanka' in data:
                rat = [data['vanka'][k]['count'] / data['star'][k]['count'] for k in data['star']
                       if data['star'][k]['status'] == 'ok' and data['vanka'][k]['status'] == 'ok' and data['star'][k]['count'] > 0]
                nfail_s = sum(1 for r in data['star'].values() if r['status'] not in ('ok', 'missing'))
                nfail_v = sum(1 for r in data['vanka'].values() if r['status'] not in ('ok', 'missing'))
                nmiss = sum(1 for r in data['vanka'].values() if r['status'] == 'missing') + sum(1 for r in data['star'].values() if r['status'] == 'missing')
                if rat:
                    lines.append(f'\nsingle column / star over the {len(rat)} runs where both converged: median {np.median(rat):.2f}, '
                                 f'mean {np.mean(rat):.2f}, min {np.min(rat):.2f}, max {np.max(rat):.2f}; not converged: star {nfail_s}, column {nfail_v}; '
                                 f'runs still missing: {nmiss}\n')
            open(os.path.join(OUT, f'summary_{mode}_{test}.md'), 'w').write('\n'.join(lines) + '\n')

            # ---- counts vs dt
            fig, ax = plt.subplots(figsize=ps.SINGLE)
            cols = ps.ramp(len(mem))
            first = 'star' if 'star' in data else 'vanka'
            for i, (lab, _, val) in enumerate(mem):
                for patch, style in [('star', ps.STAR), ('vanka', ps.COLUMN)]:
                    if patch not in data:
                        continue
                    xs = [dt for dt in S['dts'] if data[patch][(val, dt)]['status'] == 'ok']
                    ys = [data[patch][(val, dt)]['count'] for dt in xs]
                    ax.plot(xs, ys, color=cols[i], label=lab if patch == first else None, **style)
            if 'star' in data and 'vanka' in data:
                ps.patch_legend(ax)
            ax.set_xscale('log')
            ax.set_ylim(bottom=0)
            ps.finish(ax, xlabel=r'$\Delta t$ (s)', ylabel='KSP iterations' if mode == 'linear' else 'total KSP iterations', ncol=2, small=True,
                      title=(f"{mode} {test} test, $\\delta = C_1/\\Delta t$" + (f", $C_1 = {S['C1']:g}$" if test != 'dt' else '')) if TITLES else None)
            fig.savefig(os.path.join(OUT, f'counts_{mode}_{test}.png')); plt.close(fig)

            # ---- histories at the paper's dt
            # error histories at the paper's dt: one panel per patch (same conventions as the paper's figures:
            # linear histories drop the converged last point, nonlinear ones keep every Newton iterate)
            for dt in S['paper_dt'][test]:
                fig, axs = plt.subplots(1, 2, figsize=ps.WIDE, sharey=True)
                for ax, patch in zip(axs, ['star', 'vanka']):
                    if patch not in data:
                        ax.set_visible(False); continue
                    for i, (lab, _, val) in enumerate(mem):
                        r = data[patch][(val, dt)]
                        style = ps.STAR if patch == 'star' else ps.COLUMN
                        if mode == 'linear' and r['error'] is not None:
                            e = r['error'][:-1] if len(r['error']) > 1 else r['error']
                            ax.semilogy(np.arange(len(e)), e, color=cols[i], label=lab, **style)
                        elif mode == 'nonlinear' and r['snes_error'] is not None and r['snes_cum'] is not None:
                            ax.semilogy(r['snes_cum'], np.clip(r['snes_error'], 1e-16, None), color=cols[i], label=lab, **style)
                    ax.set_title('star patch' if patch == 'star' else 'single-column patch')
                    ps.finish(ax, xlabel='KSP iteration' if mode == 'linear' else 'cumulative KSP iterations', legend=(patch == 'star' or 'star' not in data), small=True)
                    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
                axs[0].set_ylabel('relative error' if mode == 'linear' else r'$\|U_k - U^*\| / \|U^*\|$')
                fig.savefig(os.path.join(OUT, f'history_{mode}_{test}_dt{dt:g}.png')); plt.close(fig)
            panels[test] = (mem, data)

        # ---- overview in the style of the paper's fig. 14
        if len(panels) == 4 and all('star' in panels[t][1] and 'vanka' in panels[t][1] for t in panels):
            fig, ax = plt.subplots(2, 2, figsize=ps.QUAD)
            (a, b), (c, d) = ax
            mem, data = panels['ar']
            cols = ps.ramp(len(mem))
            for i, (lab, _, val) in enumerate(mem):
                for patch, style in [('star', ps.STAR), ('vanka', ps.COLUMN)]:
                    xs = [dt for dt in S['dts'] if data[patch][(val, dt)]['status'] == 'ok']
                    a.plot(xs, [data[patch][(val, dt)]['count'] for dt in xs], color=cols[i], label=lab if patch == 'star' else None, **style)
            ps.patch_legend(a, column_label='single column')
            a.set_xscale('log'); a.set_title(f"(a) aspect ratio, $C_1 = {ps.sci(S['C1'])[1:-1]}$")
            a.legend(ncol=2, loc='upper right', **ps.LEGEND_SMALL)
            a.set_ylim(0, (1.35 if mode == 'nonlinear' else 1.8) * max(max(l.get_ydata()) for l in a.get_lines() if len(l.get_ydata())))   # room for the legend
            dt50 = 50.
            for test, mk, colr in [('dx', 'o', 'C0'), ('dz', 's', 'C1')]:
                mem, data = panels[test]
                ns = [val for (_, _, val) in mem]
                for patch, ls, mfc in [('star', '-', colr), ('vanka', '--', 'white')]:
                    ys = [data[patch][(n, dt50)]['count'] for n in ns]
                    b.plot(ns, ys, ls=ls, marker=mk, color=colr, mfc=mfc, label=f"{'star' if patch == 'star' else 'single column'}, {'$N_x$' if test == 'dx' else '$N_z$'}")
            b.set_title(f'(b) mesh resolution at $\\Delta t = 50$ s'); b.legend(ncol=2, **ps.LEGEND_SMALL); b.set_xticks([10, 20, 40, 60, 80])
            b.set_ylim(0, 1.35 * max(max(l.get_ydata()) for l in b.get_lines() if len(l.get_ydata())))
            mem, data = panels['dt']
            cols = ps.ramp(len(mem))
            for i, (lab, _, val) in enumerate(mem):
                for patch, style in [('star', ps.STAR), ('vanka', ps.COLUMN)]:
                    xs = [dt for dt in S['dts'] if data[patch][(val, dt)]['status'] == 'ok']
                    c.plot(xs, [data[patch][(val, dt)]['count'] for dt in xs], color=cols[i], label=lab if patch == 'star' else None, **style)
            c.set_xscale('log'); c.set_title(r'(c) shift $\delta = C_1/\Delta t$, nominal geometry'); c.legend(ncol=2, loc='upper right', **ps.LEGEND_SMALL)
            c.set_ylim(0, 1.6 * max(max(l.get_ydata()) for l in c.get_lines() if len(l.get_ydata())))   # room for the legend
            pts = []
            for test in ['dx', 'dz', 'ar']:
                mem, data = panels[test]
                for (lab, key, val) in mem:
                    for dt in S['dts']:
                        rs, rv = data['star'][(val, dt)], data['vanka'][(val, dt)]
                        if rs['status'] == 'ok' and rv['status'] == 'ok' and rs['count'] > 0:
                            nx = val if test == 'dx' else 20
                            Lz = val if test == 'ar' else LZ
                            m = 1 + (dt * N_FREQ) ** 2 / 4
                            rho1 = m * Lz ** 2 / (np.pi ** 2 * (LX / (4 * nx)) ** 2)
                            pts.append((test, rho1, rv['count'] / rs['count']))
            for test, mk, colr, name in [('dx', 'o', 'C0', 'dx sweep'), ('dz', 's', 'C1', 'dz sweep'), ('ar', '^', 'C2', 'aspect-ratio sweep')]:
                q = [p for p in pts if p[0] == test]
                d.scatter([p[1] for p in q], [p[2] for p in q], marker=mk, color=colr, label=name, s=20, edgecolor='white', linewidth=0.5)
            d.axhline(1, color=ps.GREY, lw=0.9); d.set_xscale('log'); d.set_xlabel(r'$\rho_1 = m L_z^2/(\pi^2 \Delta x^2)$'); d.set_ylabel('single column / star')
            d.set_title(f'(d) all runs where both converged ({len(pts)})'); d.legend(**ps.LEGEND_SMALL)
            for axis in (a, c):
                axis.set_xlabel(r'$\Delta t$ (s)')
            b.set_xlabel(r'$N_x$ (coarse mesh) or $N_z$')
            for axis in (a, b, c):
                axis.set_ylabel('KSP iterations' if mode == 'linear' else 'total KSP iterations')
            if TITLES:
                fig.suptitle(f"{mode} code, the paper's sweeps ($\\delta = C_1/\\Delta t$): star (solid) against single-column patch (dashed, open)")
            fig.tight_layout(h_pad=1.4, w_pad=1.2)
            fig.savefig(os.path.join(OUT, f'overview_{mode}.png')); plt.close(fig)

    with open(os.path.join(OUT, 'all_counts.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['mode', 'test', 'patch', 'member', 'dt', 'count', 'newton_steps', 'status'])
        for r in allrows:
            w.writerow([r['mode'], r['test'], r['patch'], r['member'], r['dt'], r['count'], r['newton'], r['status']])
    print('written to', OUT)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--modes', nargs='+', default=['linear', 'nonlinear'])
    ap.add_argument('--tests', nargs='+', default=['dt', 'ar', 'dx', 'dz'])
    ap.add_argument('--titles', action='store_true', help='working titles on the figures (off for the paper)')
    a = ap.parse_args()
    TITLES = a.titles
    main(a.modes, a.tests)
