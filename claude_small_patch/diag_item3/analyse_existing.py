'''
Item 3: what fig. 14(a) of the paper (claude_column_patch/paper_figure.py, panel (a)) was run with, and how the
fixed-C_1 sweeps of the paper (delta = C_1/dt) relate to the fixed-K sweeps (delta = K delta_thr(dt, eps, N)).

Reads only existing data:
    ../../data_claude_paper_nonlinear_dt_patchstar      K sweep, star, no EW   (nonlinear, the fig. 14 setting)
    ../../data_claude_paper_nonlinear_dt_ew             K sweep, star, EW on
    ../../data_claude_paper_linear_dt                   K sweep, linear code, rtol 1e-8
    ../../data_claude_paper_nonlinear_ar_patchstar      K = 1 aspect-ratio sweep (fig. 14(a), star)
    ../../data_claude_paper_nonlinear_ar_patchvanka     K = 1 aspect-ratio sweep (fig. 14(a), column)
    ../../data_nonlinear_dt_deg2_ew                     the paper's nonlinear shift sweep (fixed C_1, EW)
    ../../data_linear_dt_deg2                           the paper's linear shift sweep (fixed C_1)
Writes (next to this file): per_step_counts.md, C1_vs_K.md, fig_C1_trajectories.png, fig_per_step.png
'''
import os, re, glob, csv, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
YR2 = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(YR2, 'Boussinesq_implicit'))
import claude_paper_shift as cps  # noqa: E402

LX, LZ, N = 3.0e5, 4.0e3, 1.0e-2


def find(d, prefix, dt, key, target):
    best, bd = None, 1e9
    for f in glob.glob(f'{d}/{prefix}_dt{float(dt)!r}_{key}*.out'):
        v = float(re.search(rf'_{key}(.+)\.out$', os.path.basename(f)).group(1))
        dd = abs(np.log(v) - np.log(target)) if (v > 0 and target > 0) else abs(v - target)
        if dd < bd:
            best, bd = f, dd
    return best if bd <= 0.2 else None


def cum(path):
    return np.atleast_1d(np.loadtxt(path)).astype(int)


out = []
out.append('# Per-Newton-step KSP counts in the fixed-K sweeps (star patch, nominal geometry)\n')
out.append('delta = K delta_thr(dt); numbers are KSP iterations of Newton step 1, 2, ... (sum = total plotted in fig. 14 / the summary tables).\n')
for folder, title in [('data_claude_paper_nonlinear_dt_patchstar', 'nonlinear, no Eisenstat-Walker (the setting of fig. 14)'),
                      ('data_claude_paper_nonlinear_dt_ew', 'nonlinear, Eisenstat-Walker on (rtol0 = 1e-2)')]:
    d = f'{YR2}/{folder}'
    rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
    Ks = sorted({float(r['K']) for r in rows})
    dts = sorted({float(r['dt']) for r in rows})
    out.append(f'\n## {title}  (`{folder}`)\n')
    out.append('| dt | eps dt N | 4 eps^2 m | rho_01(K=1) | rho_10(K=1) | ' + ' | '.join(f'K={K:g}' for K in Ks) + ' |')
    out.append('|---|---|---|---|---|' + '---|' * len(Ks))
    for dt in dts:
        cells = []
        for K in Ks:
            r = [x for x in rows if float(x['K']) == K and float(x['dt']) == dt][0]
            f = find(d, 'snes_ksp_cum', dt, 'shift', float(r['delta']))
            if f is None:
                cells.append(r['status']); continue
            c = cum(f); per = np.diff(c)
            cells.append('+'.join(str(int(p)) for p in per) + (f' = {int(c[-1])}' if len(per) > 1 else '') + ('' if r['status'] == 'ok' else f' ({r["status"]})'))
        g4 = cps.group_4eps2m(dt, LX, LZ, N)
        rho10 = 1.0 / (1.0 + max(1.0, g4) / 1.0 * 2.0 / 1.0) if False else None
        # rho_10 = delta/(delta + sigma_10), sigma_10 = 4 pi^2 dt / Lx^2 ; delta = delta_thr(dt)
        sigma10 = 4 * np.pi ** 2 * dt / LX ** 2
        rho10 = cps.shift_threshold(dt, LX, LZ, N) / (cps.shift_threshold(dt, LX, LZ, N) + sigma10)
        out.append(f'| {dt:g} | {cps.group_eps_dt_N(dt, LX, LZ, N):.3g} | {g4:.3g} | {cps.rho_01(1.0, dt, LX, LZ, N):.3g} | {rho10:.3g} | ' + ' | '.join(cells) + ' |')

# linear code K sweep for comparison
d = f'{YR2}/data_claude_paper_linear_dt'
rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
Ks = sorted({float(r['K']) for r in rows}); dts = sorted({float(r['dt']) for r in rows})
out.append(f'\n## linear code LB_time_slice.py, GMRES iterations to rtol 1e-8  (`data_claude_paper_linear_dt`)\n')
out.append('| dt | ' + ' | '.join(f'K={K:g}' for K in Ks) + ' |')
out.append('|---|' + '---|' * len(Ks))
for dt in dts:
    cells = []
    for K in Ks:
        r = [x for x in rows if float(x['K']) == K and float(x['dt']) == dt][0]
        cells.append(r['ksp_its'] if r['status'] == 'ok' else r['status'])
    out.append(f'| {dt:g} | ' + ' | '.join(str(int(float(c))) if c.replace('.', '').isdigit() else c for c in cells) + ' |')

# aspect-ratio sweep of fig 14(a): per-step counts, star and column
for patch in ['star', 'vanka']:
    d = f'{YR2}/data_claude_paper_nonlinear_ar_patch{patch}'
    rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
    hs = sorted({float(r['height']) for r in rows}, reverse=True); dts = sorted({float(r['dt']) for r in rows})
    out.append(f'\n## fig. 14(a) data, {"star" if patch == "star" else "single-column"} patch, K = 1: per-Newton-step counts  (`{os.path.basename(d)}`)\n')
    out.append('| dt | ' + ' | '.join(f'Lz={h/1000:g} km' for h in hs) + ' |')
    out.append('|---|' + '---|' * len(hs))
    for dt in dts:
        cells = []
        for h in hs:
            r = [x for x in rows if float(x['height']) == h and float(x['dt']) == dt][0]
            f = find(d, 'snes_ksp_cum', dt, 'ar', h / LX)
            if f is None:
                cells.append(r['status']); continue
            c = cum(f); per = np.diff(c)
            cells.append('+'.join(str(int(p)) for p in per) + (f'={int(c[-1])}' if len(per) > 1 else ''))
        out.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')

# the paper's fixed-C1 sweeps expressed in K
out.append('\n# The paper\'s fixed-C_1 sweeps in terms of K = C_1 / (dt delta_thr(dt))\n')
out.append('K(dt) for a fixed C_1 on the nominal geometry (weak regime: K = C_1 L_x^2/(2 pi^2 dt^2); strong regime: K = 2 C_1 m L_z^2/(pi^2 dt^2)):\n')
C1s = [1e-3, 1e-4, 1e-5, 1e-6]
dts_all = [1., 5., 10., 50., 100., 500., 1000., 5000., 10000., 50000.]
out.append('| dt | regime | ' + ' | '.join(f'C_1={c:g}' for c in C1s) + ' |')
out.append('|---|---|' + '---|' * len(C1s))
for dt in dts_all:
    out.append(f'| {dt:g} | {cps.regime(dt, LX, LZ, N)} | ' + ' | '.join(f'{c/(dt*cps.shift_threshold(dt, LX, LZ, N)):.3g}' for c in C1s) + ' |')


def paper_sweep(folder, key, mode):
    d = f'{YR2}/{folder}'
    rows = {}
    for f in glob.glob(f'{d}/{"snes_ksp_cum" if mode == "nonlinear" else "error"}_dt*_{key}*.out'):
        m = re.search(rf'dt([^_]+)_{key}(.+)\.out$', os.path.basename(f))
        dt, sh = float(m.group(1)), float(m.group(2))
        a = np.atleast_1d(np.loadtxt(f))
        rows[(dt, sh)] = int(a[-1]) if mode == 'nonlinear' else len(a) - 1
    return rows


nl = paper_sweep('data_nonlinear_dt_deg2_ew', 'shift', 'nonlinear')
li = paper_sweep('data_linear_dt_deg2', 'shift', 'linear')
out.append('\n## The paper\'s nonlinear shift sweep (fig. 12 data, `data_nonlinear_dt_deg2_ew`, EW on): total KSP, with K(dt) of each run\n')
C1_nl = [1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6]
out.append('| dt | ' + ' | '.join(f'C_1={c:g}' for c in C1_nl) + ' |')
out.append('|---|' + '---|' * len(C1_nl))
for dt in dts_all:
    cells = []
    for c in C1_nl:
        sh = c / dt
        hit = [v for (ddt, s), v in nl.items() if abs(ddt - dt) < 1e-9 and abs(np.log(s) - np.log(sh)) < 0.2]
        K = c / (dt * cps.shift_threshold(dt, LX, LZ, N))
        cells.append((f'{hit[0]}' if hit else 'missing') + f' (K={K:.2g})')
    out.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')
out.append('\n## The paper\'s linear shift sweep (fig. 7 data, `data_linear_dt_deg2`, rtol 1e-8): KSP iterations, with K(dt)\n')
C1_li = [1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-5]
dts_li = [1., 5., 10., 50., 100., 500., 1000., 5000., 10000.]
out.append('| dt | ' + ' | '.join(f'C_1={c:g}' for c in C1_li) + ' |')
out.append('|---|' + '---|' * len(C1_li))
for dt in dts_li:
    cells = []
    for c in C1_li:
        sh = c / dt
        hit = [v for (ddt, s), v in li.items() if abs(ddt - dt) < 1e-9 and abs(np.log(s) - np.log(sh)) < 0.2]
        K = c / (dt * cps.shift_threshold(dt, LX, LZ, N))
        cells.append((f'{hit[0]}' if hit else 'missing') + f' (K={K:.2g})')
    out.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')
open(os.path.join(HERE, 'per_step_counts_and_C1_vs_K.md'), 'w').write('\n'.join(out) + '\n')

# ---- figure: iteration map of the K sweep with the fixed-C1 trajectories of the paper overlaid
fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
for ax, (folder, title) in zip(axs, [('data_claude_paper_nonlinear_dt_patchstar', 'nonlinear code, star patch, no EW (fig. 14 setting)'),
                                      ('data_claude_paper_linear_dt', 'linear code (rtol 1e-8)')]):
    d = f'{YR2}/{folder}'
    rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
    Ks = sorted({float(r['K']) for r in rows}); dts = sorted({float(r['dt']) for r in rows})
    Z = np.full((len(Ks), len(dts)), np.nan)
    for r in rows:
        if r['status'] == 'ok':
            v = r['total_ksp'] if 'nonlinear' in folder else r['ksp_its']
            Z[Ks.index(float(r['K'])), dts.index(float(r['dt']))] = float(v)
    X, Y = np.meshgrid(dts, Ks)
    pm = ax.pcolormesh(X, Y, Z, shading='nearest', cmap='viridis', norm=matplotlib.colors.LogNorm(vmin=3, vmax=150))
    for i, K in enumerate(Ks):
        for j, dt in enumerate(dts):
            ax.text(dt, K, f'{int(Z[i, j])}' if np.isfinite(Z[i, j]) else 'x', ha='center', va='center', fontsize=7,
                    color='white' if (np.isfinite(Z[i, j]) and Z[i, j] < 40) else 'black')
    dtf = np.logspace(0, np.log10(max(dts)), 200)
    for c, ls in [(1e-3, '-'), (1e-5, '--'), (1e-6, ':')] if 'nonlinear' in folder else [(1e-3, '-'), (1e-4, '--'), (1e-5, ':')]:
        ax.plot(dtf, [c / (t * cps.shift_threshold(t, LX, LZ, N)) for t in dtf], 'r' + ls, lw=1.5, label=rf'paper: $C_1 = {c:g}$')
    ax.axvline(cps.crossover_dt(LX, LZ, N), color='w', ls='--', lw=0.8)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel(r'$\Delta t$ (s)'); ax.set_ylabel(r'$K = \delta/\delta_{thr}$')
    ax.set_ylim(min(Ks) / 3, max(Ks) * 3); ax.set_title(title, fontsize=10); ax.legend(fontsize=8, loc='lower left')
fig.colorbar(pm, ax=axs, label='KSP iterations (total over Newton steps)', shrink=0.9)
fig.suptitle(r'Fixed-$K$ sweeps (cells) and the paper\'s fixed-$C_1$ sweeps (red lines: $K(\Delta t) = C_1/(\Delta t\,\delta_{thr})$)', fontsize=11)
fig.savefig(os.path.join(HERE, 'fig_C1_trajectories.png'), dpi=160, bbox_inches='tight')

# ---- figure: per-step counts vs dt, K = 0.1, 1, 10 (star, no EW) and the linear code
fig, ax = plt.subplots(figsize=(7.5, 4.8))
d = f'{YR2}/data_claude_paper_nonlinear_dt_patchstar'
rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
for K, col in [(0.1, 'C0'), (1.0, 'C1'), (10.0, 'C3')]:
    dts = sorted({float(r['dt']) for r in rows}); s1, s2, tot = [], [], []
    for dt in dts:
        r = [x for x in rows if float(x['K']) == K and float(x['dt']) == dt][0]
        f = find(d, 'snes_ksp_cum', dt, 'shift', float(r['delta']))
        c = cum(f); per = np.diff(c)
        s1.append(per[0]); s2.append(per[1] if len(per) > 1 else np.nan); tot.append(c[-1])
    ax.plot(dts, s1, '-o', color=col, label=f'K={K:g}: Newton step 1')
    ax.plot(dts, s2, '--s', color=col, mfc='white', label=f'K={K:g}: Newton step 2')
d = f'{YR2}/data_claude_paper_linear_dt'
rows = list(csv.DictReader(open(f'{d}/summary_iterations.csv')))
dts = sorted({float(r['dt']) for r in rows})
ax.plot(dts, [float([x for x in rows if float(x['K']) == 1.0 and float(x['dt']) == dt][0]['ksp_its']) for dt in dts], '-^', color='k', label='linear code, K=1 (rtol 1e-8)')
ax.axvline(cps.crossover_dt(LX, LZ, N), color='grey', ls='--', lw=0.8)
ax.set_xscale('log'); ax.set_xlabel(r'$\Delta t$ (s)'); ax.set_ylabel('KSP iterations per Newton step'); ax.grid(alpha=0.3)
ax.set_title('Nonlinear code, star patch, nominal geometry, no EW: iterations per Newton step', fontsize=10); ax.legend(fontsize=8)
fig.savefig(os.path.join(HERE, 'fig_per_step.png'), dpi=160, bbox_inches='tight')
print('written', os.path.join(HERE, 'per_step_counts_and_C1_vs_K.md'))
