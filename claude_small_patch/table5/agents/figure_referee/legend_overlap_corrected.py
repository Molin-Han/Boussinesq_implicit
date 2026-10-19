'''Which plotted markers / line segments of the restyled figures lie under the legend box?  Re-runs the restyle
functions with Figure.savefig patched to (i) report overlaps, (ii) write into this folder only.'''
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch')
import restyle_paper_figures as R
import paper_style as ps
OUT = os.path.dirname(os.path.abspath(__file__)) + '/overlap_corrected'
os.makedirs(OUT, exist_ok=True)
_orig = Figure.savefig

def patched(self, fname, *a, **k):
    self.canvas.draw()
    r = self.canvas.get_renderer()
    for ax in self.axes:
        leg = ax.get_legend()
        if leg is None:
            continue
        bb = leg.get_window_extent(r)
        hits = []
        for ln in ax.get_lines():
            lab = ln.get_label()
            if lab.startswith('_'):
                continue
            xy = np.column_stack(ln.get_data())
            xy = xy[np.all(np.isfinite(xy), axis=1)]
            if len(xy) == 0:
                continue
            P = ax.transData.transform(xy)
            ms = ln.get_markersize() * self.dpi / 72 / 2
            inside = [(x, y) for (x, y), (px, py) in zip(xy, P)
                      if bb.x0 - ms <= px <= bb.x1 + ms and bb.y0 - ms <= py <= bb.y1 + ms]
            # segments crossing the box (sampled)
            seg = 0
            for (p, q) in zip(P[:-1], P[1:]):
                t = np.linspace(0, 1, 200)[:, None]
                S = p + t * (q - p)
                seg += np.sum((S[:, 0] >= bb.x0) & (S[:, 0] <= bb.x1) & (S[:, 1] >= bb.y0) & (S[:, 1] <= bb.y1))
            if inside or seg:
                hits.append(f'{lab}: markers under legend {[(float(f"{x:.4g}"), float(f"{y:.3g}")) for x, y in inside]}, '
                            f'segment samples under legend {seg}')
        print(f'{os.path.basename(str(fname))}: legend frame alpha {leg.get_frame().get_alpha()}; '
              + ('; '.join(hits) if hits else 'no overlap'))
    return _orig(self, os.path.join(OUT, os.path.basename(str(fname))), *a, **k)

Figure.savefig = patched
sys.argv = ['restyle_paper_figures.py', '--set', 'paper', '--out', OUT, '--corrected']
R.main()
