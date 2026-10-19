'''
Figure for the paper section on the single-column patch (Patch-Preconditioner-for-Boussinesq).

    python paper_figure.py            writes ../../Patch-Preconditioner-for-Boussinesq/Image_patch/column_patch_sweeps.png

Panels, all from the sweep folders ../../data_claude_paper_nonlinear_<test>_patch{star,vanka}:
  (a) aspect-ratio sweep, total FGMRES iterations against dt, star (solid) and single column (dashed);
  (b) mesh sweeps at dt = 50 s: iterations against N_x (dx sweep) and N_z (dz sweep);
  (c) shift sweep on the nominal geometry, K in {0.1, 1, 10, 100, 1000} (runs that did not converge are omitted);
  (d) ratio single column / star of all K = 1 runs against rho_1 = m L_z^2 / (pi^2 dx^2).
'''
import os
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..'))
OUT = os.path.abspath(os.path.join(HERE, '..', '..', 'Patch-Preconditioner-for-Boussinesq', 'Image_patch'))
N_FREQ = 1.0e-2

plt.rcParams.update({'font.size': 10, 'axes.titlesize': 10.5, 'legend.fontsize': 7.5, 'axes.labelsize': 10,
                     'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5, 'lines.linewidth': 1.6, 'lines.markersize': 4})
RAMP = ['#86b6ef', '#5598e7', '#2a78d6', '#256abf', '#1c5cab', '#104281', '#0d366b']   # one hue, light -> dark
CAT = {'dx': '#2a78d6', 'dz': '#eb6834', 'ar': '#1baf7a'}
GREY = '#6b6a63'


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


fig, ax = plt.subplots(2, 2, figsize=(7.4, 5.9))
(a, b), (c, d) = ax

# (a) aspect ratio
S, V = table('ar', 'star'), table('ar', 'vanka')
dts = sorted({k[0] for k in S})
heights = sorted({k[1] for k in S}, reverse=True)
for i, h in enumerate(heights):
    col = RAMP[round(i * (len(RAMP) - 1) / (len(heights) - 1))]
    lab = f'{h / 1000:g} km'
    xs, ys = series(S, h, dts); a.plot(xs, ys, '-', marker='o', color=col, label=lab if i in (0, 3, 6) else None)
    xs, ys = series(V, h, dts); a.plot(xs, ys, '--', marker='o', mfc='white', color=col)
a.set_xscale('log'); a.set_xlabel(r'$\Delta t$ (s)'); a.set_ylabel('total FGMRES iterations')
a.set_title(r'(a) aspect ratio, $L_z$ from 12 km (light) to 0.25 km (dark)')
a.set_ylim(0, 50)
h1, = a.plot([], [], '-', marker='o', color=GREY, label='star patch')
h2, = a.plot([], [], '--', marker='o', mfc='white', color=GREY, label='single column')
a.legend(loc='upper left', ncol=2, handlelength=2.2)

# (b) mesh sweeps at dt = 50
for test, key, mk in (('dx', 'N_x', 'o'), ('dz', 'N_z', 's')):
    S, V = table(test, 'star'), table(test, 'vanka')
    ns = sorted({k[1] for k in S})
    b.plot(ns, [S[(50.0, n)][0] for n in ns], '-', marker=mk, color=CAT[test], label=f'star, ${key}$')
    b.plot(ns, [V[(50.0, n)][0] for n in ns], '--', marker=mk, mfc='white', color=CAT[test], label=f'single column, ${key}$')
b.set_xlabel(r'$N_x$ (coarse mesh) or $N_z$'); b.set_ylabel('total FGMRES iterations')
b.set_title(r'(b) mesh resolution at $\Delta t = 50$ s'); b.set_ylim(0, 50); b.set_xticks([10, 20, 40, 60, 80])
b.legend(loc='upper left', ncol=2)

# (c) shift sweep
S, V = table('dt', 'star'), table('dt', 'vanka')
Ks = sorted({k[1] for k in S})
for i, K in enumerate(Ks):
    col = RAMP[round(i * (len(RAMP) - 1) / (len(Ks) - 1))]
    xs, ys = series(S, K, dts); c.plot(xs, ys, '-', marker='o', color=col, label=f'$K = {K:g}$')
    xs, ys = series(V, K, dts); c.plot(xs, ys, '--', marker='o', mfc='white', color=col)
c.set_xscale('log'); c.set_xlabel(r'$\Delta t$ (s)'); c.set_ylabel('total FGMRES iterations')
c.set_title(r'(c) shift $\delta = K\,\delta_{\mathrm{thr}}$, nominal geometry'); c.set_ylim(0, 140)
c.legend(loc='upper left', ncol=3)

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
    d.scatter([p[1] for p in q], [p[2] for p in q], s=22, marker=mk, color=CAT[test], edgecolor='white', linewidth=0.5, label=name)
d.axhline(1.0, color=GREY, lw=1.0)
d.set_xscale('log'); d.set_xlabel(r'$\rho_1 = m L_z^2 / (\pi^2 \Delta x^2)$'); d.set_ylabel('iterations, single column / star')
d.set_title(r'(d) all $K = 1$ runs against $\rho_1$'); d.set_ylim(0, 11.5)
d.legend(loc='upper right')

for axis in ax.flat:
    axis.grid(alpha=0.3)
fig.tight_layout(h_pad=1.4, w_pad=1.2)
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, 'column_patch_sweeps.png'), dpi=220)
print('written', os.path.join(OUT, 'column_patch_sweeps.png'), len(pts), 'points in (d)')
