'''Trial legend placements for the three restyled figures whose legend hides data. Re-runs the restyle functions;
just before saving, the legend is rebuilt with the trial settings, then the overlap test of legend_overlap.py runs.
Writes only into ./fix_trial/.'''
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import restyle_paper_figures as R
import paper_style as ps
OUT = os.path.dirname(os.path.abspath(__file__)) + '/fix_trial'
os.makedirs(OUT, exist_ok=True)
_orig = Figure.savefig
TRIAL = {}

def overlap_report(fig, name):
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    ax = fig.axes[0]; leg = ax.get_legend(); bb = leg.get_window_extent(r)
    hits = []
    for ln in ax.get_lines():
        if ln.get_label().startswith('_'):
            continue
        xy = np.column_stack(ln.get_data()); xy = xy[np.all(np.isfinite(xy), axis=1)]
        if not len(xy):
            continue
        P = ax.transData.transform(xy); ms = ln.get_markersize() * fig.dpi / 72 / 2
        ins = [tuple(map(float, p)) for p, q in zip(xy, P) if bb.x0 - ms <= q[0] <= bb.x1 + ms and bb.y0 - ms <= q[1] <= bb.y1 + ms]
        seg = 0
        for p, q in zip(P[:-1], P[1:]):
            S = p + np.linspace(0, 1, 200)[:, None] * (q - p)
            seg += np.sum((S[:, 0] >= bb.x0) & (S[:, 0] <= bb.x1) & (S[:, 1] >= bb.y0) & (S[:, 1] <= bb.y1))
        if ins or seg:
            hits.append(f'{ln.get_label()}: markers {ins}, segment samples {seg}')
    return 'no overlap' if not hits else '; '.join(hits)

def patched(self, fname, *a, **k):
    name = os.path.basename(str(fname))
    for tag, fix in TRIAL.get(name, []):
        fix(self)
        out = os.path.join(OUT, name.replace('.png', f'__{tag}.png'))
        print(f'{name} [{tag}]: {overlap_report(self, name)}')
        _orig(self, out, *a, **k)
    return None

def relabel_short_ar(fig):
    ax = fig.axes[0]
    for ln in ax.get_lines():
        lab = ln.get_label()
        if lab.startswith('$L_z'):
            ln.set_label(lab.split(',')[0])

def legend(**kw):
    def f(fig):
        ax = fig.axes[0]
        if ax.get_legend():
            ax.get_legend().remove()
        ax.legend(**kw)
    return f

def ylim_top(top):
    def f(fig):
        fig.axes[0].set_ylim(top=top)
    return f

def chain(*fs):
    def f(fig):
        for g in fs:
            g(fig)
    return f

S = ps.LEGEND_SMALL
TRIAL['error_AR_Robust_dt50.0.png'] = [
    ('shortlabels', chain(relabel_short_ar, legend(loc='upper right'))),
]
TRIAL['error_shift.png'] = [
    ('ncol2_ylim120', chain(legend(loc='upper right', ncol=2), ylim_top(120))),
    ('ncol3_ylim115', chain(legend(loc='upper center', ncol=3), ylim_top(115))),
]
TRIAL['snes_error_shift_dt500.0.png'] = [
    ('lowerleft_ncol1', legend(loc='lower left', fontsize=6.5, handlelength=1.6, borderpad=0.3, labelspacing=0.25)),
    ('below_ncol3', legend(loc='upper center', bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize=6.5, handlelength=1.6,
                           borderpad=0.3, labelspacing=0.25, columnspacing=0.8)),
    ('outside_right', legend(loc='upper left', bbox_to_anchor=(1.01, 1.0), fontsize=6.5, handlelength=1.6,
                             borderpad=0.3, labelspacing=0.25)),
]
Figure.savefig = patched
ps.apply()
R.linear_history(R.folder('paper', 'linear', 'ar'), 'ar', 50.0, False, OUT, 'error_AR_Robust_dt50.0.png', ps.SINGLE)
R.linear_shift_sweep(R.folder('paper', 'linear', 'dt'), False, OUT, 'error_shift.png', ps.SINGLE)
R.nonlinear_history(R.folder('paper', 'nonlinear', 'dt'), 'dt', 500.0, False, OUT, 'snes_error_shift_dt500.0.png', ps.HALF)
