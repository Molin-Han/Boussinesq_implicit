'''
The LSR report's iteration figures from the 2026-11-16 re-runs, drawn with the paper's figure code
(make_paper_figures_gmres6.py: same data, counting convention and style).

The LSR includes the same figures as the paper under its own folder names, plus two the paper does not have:
  Image_new/error_{dx,dz,AR}_Robust_dt50.0.png, error_shift.png    linear; placed alone at 0.8-0.85 \\linewidth,
                                                                   so drawn at the paper's single-figure size
  Image_nonlinear/snes_error_{dx,dz}_dt50.0.png, snes_error_AR_dt100.0.png, error_shift_nonlinear.png
  Image_nonlinear/snes_error_AR_dt5000.0.png                       aspect-ratio histories at dt = 5000 s
  Image_nonlinear/snes_error_shift_dt500.0.png                     shift-sweep histories at dt = 500 s
(nonlinear ones sit in half-width minipages, so they keep the paper's half size).

    python make_lsr_figures_gmres6.py [--data DIR] [--lsr DIR] [--out DIR]
'''
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_paper_figures_gmres6 as m  # noqa: E402
import paper_style as ps  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def nonlinear_shift_history(data, dt, out):
    '''histories of the nonlinear shift sweep at one time step (the runs passed the stage shift C_1/(2 dt))'''
    fig, ax = plt.subplots(figsize=ps.HALF)
    mem = m.members('nonlinear', 'dt')
    d = m.folder(data, 'nonlinear', 'dt', 'star')
    for (lab, key, c1, _), col in zip(mem, ps.ramp(len(mem))):
        tgt = m.shift_target('nonlinear', c1, dt)
        se = m.load(d, 'snes_error', dt, key, tgt)
        cu = m.load(d, 'snes_ksp_cum', dt, key, tgt)
        if se is None or cu is None or se.size == 0:
            continue
        ax.semilogy(cu, np.clip(se, m.FLOOR, None), color=col, label=lab, **ps.STAR)
    ps.finish(ax, xlabel='cumulative FGMRES iterations', ylabel=r'$\|U_k - U^*\| / \|U^*\|$', legend=False)
    ax.set_ylim(top=1e4)   # room above the curves for the seven-entry legend
    ax.legend(loc='upper right', ncol=2, **dict(m.LEG8, fontsize=7, columnspacing=0.6, handlelength=1.2))
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.savefig(out); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default='/Users/molin_han/Code_Repo/data/paper_robustness_gmres6')
    ap.add_argument('--lsr', default=os.path.abspath(os.path.join(HERE, '..', '..', '..', 'Tex_Repo', 'LSR_Copied_from_ESA')))
    ap.add_argument('--out', default=None, help='write all figures into this one folder instead of the LSR tree')
    a = ap.parse_args()
    path = (lambda sub, name: os.path.join(a.out, name)) if a.out else (lambda sub, name: os.path.join(a.lsr, sub, name))
    if a.out:
        os.makedirs(a.out, exist_ok=True)
    ps.apply()
    recs = m.all_records(a.data)
    half = ps.HALF
    ps.HALF = ps.SINGLE      # the linear figures stand alone in the LSR
    try:
        for test, name in (('dx', 'error_dx_Robust_dt50.0.png'), ('dz', 'error_dz_Robust_dt50.0.png'),
                           ('ar', 'error_AR_Robust_dt50.0.png')):
            m.linear_history(a.data, test, path('Image_new', name))
        m.linear_shift_sweep(recs, path('Image_new', 'error_shift.png'))
    finally:
        ps.HALF = half
    m.nonlinear_history(a.data, 'dx', 50.0, path('Image_nonlinear', 'snes_error_dx_dt50.0.png'))
    m.nonlinear_history(a.data, 'dz', 50.0, path('Image_nonlinear', 'snes_error_dz_dt50.0.png'))
    m.nonlinear_history(a.data, 'ar', 100.0, path('Image_nonlinear', 'snes_error_AR_dt100.0.png'))
    m.nonlinear_history(a.data, 'ar', 5000.0, path('Image_nonlinear', 'snes_error_AR_dt5000.0.png'))
    nonlinear_shift_history(a.data, 500.0, path('Image_nonlinear', 'snes_error_shift_dt500.0.png'))
    m.nonlinear_shift_sweep(recs, path('Image_nonlinear', 'error_shift_nonlinear.png'))
    print('LSR figures written to', a.out or a.lsr)


if __name__ == '__main__':
    main()
