'''
Fig. 14 of the paper (column_patch_sweeps.png) regenerated in the shared style (paper_style.py); a copy of
claude_column_patch/paper_figure.py with the colours, markers and fonts replaced.  Same data, same panels:

    python restyle_fig14.py [--out DIR]     default DIR = ../../Patch-Preconditioner-for-Boussinesq/Image_restyled

  (a) aspect-ratio sweep, total FGMRES iterations against dt, star (solid, filled) and single column (dashed, open);
  (b) mesh sweeps at dt = 50 s: iterations against N_x (dx sweep) and N_z (dz sweep);
  (c) shift sweep on the nominal geometry, K in {0.1, 1, 10, 100, 1000} (runs that did not converge are omitted);
  (d) ratio single column / star of all K = 1 runs against rho_1 = m L_z^2 / (pi^2 dx^2).
'''
import os
import sys
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_style as ps  # noqa: E402
DATA = os.path.abspath(os.path.join(HERE, '..', '..'))
ap = argparse.ArgumentParser()
ap.add_argument('--out', default=os.path.join(DATA, 'Patch-Preconditioner-for-Boussinesq', 'Image_restyled'))
OUT = ap.parse_args().out
N_FREQ = 1.0e-2
ps.apply()

CAT = {'dx': 'C0', 'dz': 'C1', 'ar': 'C2'}
GREY = ps.GREY


def read(test, patch):
    path = os.path.join(DATA, f'data_claude_paper_nonlinear_{test}_patch{patch}', 'summary_iterations.csv')
    with open(path) as f:
        return list(csv.DictReader(f))


def member(test, r):
    return {'dx': int(r['nx']), 'dz': int(r['nz']), 'ar': float(r['height']), 'dt': float(r['K'])}[test]


def table(test, patch):
    out = {}
    for r in read(test, patch):
        tot = r['total_ksp']
        out[(float(r['dt']), member(test, r))] = (int(float(tot)) if tot not in ('', 'nan') else None, r['status'])
    return out


def series(T, m, dts):
    xs = [dt for dt in dts if T.get((dt, m)) and T[(dt, m)][1] == 'ok' and T[(dt, m)][0]]
    return xs, [T[(dt, m)][0] for dt in xs]


fig, ax = plt.subplots(2, 2, figsize=ps.QUAD)
(a, b), (c, d) = ax

# (a) aspect ratio
S, V = table('ar', 'star'), table('ar', 'vanka')
dts = sorted({k[0] for k in S})
heights = sorted({k[1] for k in S}, reverse=True)
cols = ps.ramp(len(heights))
for i, h in enumerate(heights):
    col = cols[i]
    lab = f'{h / 1000:g} km'
    xs, ys = series(S, h, dts); a.plot(xs, ys, color=col, label=lab if i in (0, 3, 6) else None, **ps.STAR)
    xs, ys = series(V, h, dts); a.plot(xs, ys, color=col, **ps.COLUMN)
a.set_xscale('log'); a.set_xlabel(r'$\Delta t$ (s)'); a.set_ylabel('total FGMRES iterations')
a.set_title(r'(a) aspect ratio, $K = 1$')
a.set_ylim(0, 62)   # room for the legend above the data (max 44)
ps.patch_legend(a, star_label='star', column_label='single column')
a.legend(loc='upper left', ncol=2, **ps.LEGEND_SMALL)

# (b) mesh sweeps at dt = 50
for test, key, mk in (('dx', 'N_x', 'o'), ('dz', 'N_z', 's')):
    S, V = table(test, 'star'), table(test, 'vanka')
    ns = sorted({k[1] for k in S})
    b.plot(ns, [S[(50.0, n)][0] for n in ns], ls='-', marker=mk, color=CAT[test], label=f'star, ${key}$')
    b.plot(ns, [V[(50.0, n)][0] for n in ns], ls='--', marker=mk, mfc='white', color=CAT[test], label=f'single column, ${key}$')
b.set_xlabel(r'$N_x$ (coarse mesh) or $N_z$'); b.set_ylabel('total FGMRES iterations')
b.set_title(r'(b) mesh resolution, $\Delta t = 50$ s, $K = 1$'); b.set_ylim(0, 58); b.set_xticks([10, 20, 40, 60, 80])
b.legend(loc='upper left', ncol=2, **ps.LEGEND_SMALL)

# (c) shift sweep
S, V = table('dt', 'star'), table('dt', 'vanka')
Ks = sorted({k[1] for k in S})
cols = ps.ramp(len(Ks))
for i, K in enumerate(Ks):
    col = cols[i]
    xs, ys = series(S, K, dts); c.plot(xs, ys, color=col, label=f'$K = {K:g}$', **ps.STAR)
    xs, ys = series(V, K, dts); c.plot(xs, ys, color=col, **ps.COLUMN)
c.set_xscale('log'); c.set_xlabel(r'$\Delta t$ (s)'); c.set_ylabel('total FGMRES iterations')
c.set_title(r'(c) shift $\delta = K\,\delta_{\mathrm{thr}}$, nominal geometry'); c.set_ylim(0, 178)   # room for the legend above the K = 100 peak (131)
c.legend(loc='upper left', ncol=3, **ps.LEGEND_SMALL)

# (d) collapse
pts = []
for test in ('dx', 'dz', 'ar'):
    S, V = table(test, 'star'), table(test, 'vanka')
    for r in read(test, 'star'):
        k = (float(r['dt']), member(test, r))
        if k in V and S[k][1] == 'ok' and V[k][1] == 'ok' and S[k][0]:
            dt, Lx, Lz, nx = float(r['dt']), float(r['length']), float(r['height']), int(r['nx'])
            m = 1.0 + (dt * N_FREQ) ** 2 / 4.0
            pts.append((test, m * Lz ** 2 / (np.pi ** 2 * (Lx / (4 * nx)) ** 2), V[k][0] / S[k][0]))
for test, mk, name in (('dx', 'o', 'dx sweep'), ('dz', 's', 'dz sweep'), ('ar', '^', 'aspect-ratio sweep')):
    q = [p for p in pts if p[0] == test]
    d.scatter([p[1] for p in q], [p[2] for p in q], s=20, marker=mk, color=CAT[test], edgecolor='white', linewidth=0.5, label=name)
d.axhline(1.0, color=GREY, lw=1.0)
d.set_xscale('log'); d.set_xlabel(r'$\rho_1 = m L_z^2 / (\pi^2 \Delta x^2)$'); d.set_ylabel('single column / star')
d.set_title(r'(d) all $K = 1$ runs against $\rho_1$'); d.set_ylim(0, 11.5)
d.legend(loc='upper right', **ps.LEGEND_SMALL)

fig.tight_layout(h_pad=1.4, w_pad=1.2)
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'column_patch_sweeps.png'))
print('written', os.path.join(OUT, 'column_patch_sweeps.png'), len(pts), 'points in (d)')
