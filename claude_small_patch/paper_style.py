'''
One matplotlib style for every figure of the paper and of the patch studies.

    import paper_style as ps
    ps.apply()                       # rcParams
    cols = ps.ramp(len(members))     # ordered colours (viridis, light end cut off so it prints on white)
    ax.plot(x, y, **ps.STAR, color=cols[i], label=...)      # star patch: solid, filled circles
    ax.plot(x, y, **ps.COLUMN, color=cols[i], label=...)    # single-column patch: dashed, open squares
    ps.finish(ax, xlabel=..., ylabel=...)                   # labels, legend, grid
    fig.savefig(path)                                       # bbox_inches='tight' and dpi come from rcParams

Figure sizes are chosen so that the fonts print close to their nominal size in the SIAM one-column layout
(text width 5.125 in): ps.SINGLE = (4.6, 3.1) for a figure placed at 0.8-0.85 linewidth (prints at 4.1-4.4 in),
ps.HALF = (2.8, 2.4) for the two-panel minipage figures of the nonlinear section (each 0.48 linewidth = 2.46 in),
ps.WIDE = (5.8, 2.9) for two panels side by side in one file, ps.QUAD = (6.4, 5.2) for four panels (linewidth).
'''
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RC = {
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 10,
    'legend.fontsize': 8,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'lines.linewidth': 1.5,
    'lines.markersize': 4,
    'lines.markeredgewidth': 1.0,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'grid.linewidth': 0.6,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.framealpha': 0.92,
    'legend.edgecolor': '0.85',
    'legend.handlelength': 2.4,
    'figure.dpi': 100,
    'savefig.dpi': 220,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.03,
    'mathtext.default': 'regular',
}
SINGLE = (4.6, 3.1)
HALF = (2.8, 2.4)
WIDE = (5.8, 2.9)
QUAD = (6.4, 5.2)
LEGEND_SMALL = dict(fontsize=7, handlelength=1.8, borderpad=0.3, labelspacing=0.3, columnspacing=1.0)
STAR = dict(ls='-', marker='o')
COLUMN = dict(ls='--', marker='s', mfc='white')
GREY = '#6b6a63'


def apply():
    plt.rcParams.update(RC)


def ramp(n, cmap='viridis', lo=0.0, hi=0.82):
    '''n ordered colours from a sequential map; the brightest end is cut off (hi < 1) for legibility on white.'''
    cm = plt.get_cmap(cmap)
    if n == 1:
        return [cm(lo)]
    return [cm(lo + (hi - lo) * i / (n - 1)) for i in range(n)]


def sci(x, digits=1):
    '''$10^{-3}$, $5\\times10^{-4}$, $1.5\\times10^{2}$ ... for legends.'''
    if x == 0:
        return '$0$'
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10.0 ** e
    if abs(m - 1.0) < 1e-9:
        return rf'$10^{{{e}}}$'
    if abs(m - round(m)) < 1e-9:
        return rf'${int(round(m))}\times10^{{{e}}}$'
    return rf'${m:.{digits}f}\times10^{{{e}}}$'


def finish(ax, xlabel=None, ylabel=None, legend=True, loc='best', ncol=1, title=None, small=False):
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if legend:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(loc=loc, ncol=ncol, **(LEGEND_SMALL if small else {}))
    ax.grid(True, which='major', alpha=RC['grid.alpha'])


def patch_legend(ax, loc='upper right', star_label='star patch', column_label='single-column patch'):
    '''Grey proxy entries for the two line styles; use together with member labels on the star curves.'''
    ax.plot([], [], color=GREY, label=star_label, **STAR)
    ax.plot([], [], color=GREY, label=column_label, **COLUMN)
