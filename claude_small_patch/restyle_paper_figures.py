'''
The ten iteration figures of the paper, regenerated in the shared style (paper_style.py), with the same file
names as in the paper so that they are drop-in replacements:

    Image_new/error_dx_Robust_dt50.0.png     linear, dt = 50 s, one curve per N_x      (fig. 4)
    Image_new/error_dz_Robust_dt50.0.png     linear, dt = 50 s, one curve per N_z      (fig. 5)
    Image_new/error_AR_Robust_dt50.0.png     linear, dt = 50 s, one curve per L_z      (fig. 6)
    Image_new/error_shift.png                linear, KSP iterations against dt per C_1 (fig. 7)
    Image_nonlinear/snes_error_dx_dt50.0.png, snes_error_dz_dt50.0.png                 (figs. 8, 9)
    Image_nonlinear/snes_error_AR_dt100.0.png, snes_error_AR_dt5000.0.png              (fig. 10)
    Image_nonlinear/snes_error_shift_dt500.0.png, error_shift_nonlinear.png            (figs. 11, 12)

Data sets (--set):
    paper   the folders the paper's figures were made from: ../data_linear_<test>_deg2, ../data_nonlinear_<test>_deg2_ew
    star    the 2026-10-16 star re-run of the same sweeps:   ../data_linear_<test>_deg2_star, ../data_nonlinear_<test>_deg2_ew_star
    vanka   the single-column patch:                          ../data_linear_<test>_deg2_vanka, ../data_nonlinear_<test>_deg2_ew_vanka

Conventions.  By default the figures show exactly what the paper's figures show (so that only the style changes):
  * linear histories drop the last monitored point (the converged state), as plot_dx.py / plot_dz.py / plot_ar.py do;
  * the linear shift sweep plots len(error) = n + 1 for a run with n KSP iterations (plot_shift.py), NaN at the cap;
  * the nonlinear histories show every Newton iterate including the converged one, as four of the paper's five
    figures do; the paper's snes_error_AR_dt5000.0 was made with the newer plot_ar.py that drops the converged
    point, and the drop-in copy does the same (drop_last per figure);
  * the nonlinear shift sweep plots len(error) = total KSP iterations + number of Newton steps (plot_shift.py).
With --corrected the counts are the iteration counts themselves (n, and the total KSP count over the Newton steps),
the linear histories keep their final point, the nonlinear histories all drop the converged point (one convention,
the one of the current plot_*.py scripts), and the output goes to <out>/corrected/.

    python restyle_paper_figures.py --set paper                 # -> ../../Patch-Preconditioner-for-Boussinesq/Image_restyled/
    python restyle_paper_figures.py --set vanka --out DIR
'''
import os
import re
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
YR2 = os.path.abspath(os.path.join(HERE, '..', '..'))
import sys
sys.path.insert(0, HERE)
import paper_style as ps  # noqa: E402

LX = 3.0e5
MAXIT = 100
LIN = dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000.], C1=1e-3, C1_list=[1., 0.1, 0.01, 1e-3, 1e-4, 1e-5],
           heights=[12000, 8000, 4000, 2000, 1000, 500, 250], nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80])
NL = dict(dts=[1., 5., 10., 50., 100., 500., 1000., 5000., 10000., 50000.], C1=1e-5,
          C1_list=[1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6], heights=[12000, 8000, 4000, 2000, 1000],
          nxs=[10, 20, 40, 60, 80], nzs=[10, 20, 40, 60, 80])
SUFFIX = {'paper': '', 'star': '_star', 'vanka': '_vanka'}


def folder(dataset, mode, test):
    return f'{YR2}/data_{mode}_{test}_deg2' + ('_ew' if mode == 'nonlinear' else '') + SUFFIX[dataset]


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


def load(d, prefix, dt, key, tgt):
    f = find(d, prefix, dt, key, tgt)
    return None if f is None else np.atleast_1d(np.loadtxt(f))


def members(mode, test, short=False):
    '''(legend label, file key, value); short labels for the half-width figures.'''
    S = LIN if mode == 'linear' else NL
    if test == 'dx':
        return [(rf'$N_x = {n}$', 'nx', n) for n in S['nxs']]
    if test == 'dz':
        return [(rf'$N_z = {n}$', 'nz', n) for n in S['nzs']]
    if test == 'ar':
        if short:
            return [(rf'$L_z = {h/1000:g}$ km', 'ar', h / LX) for h in S['heights']]
        return [(rf'$L_z = {h/1000:g}$ km, $\varepsilon = {ps.sci(h / LX, 1)[1:-1]}$', 'ar', h / LX) for h in S['heights']]
    return [(rf'$C_1 = {ps.sci(c)[1:-1]}$', 'shift', c) for c in S['C1_list']]


# ---------------------------------------------------------------------------------------------------------------
def linear_history(d, test, dt, corrected, out, name, size):
    fig, ax = plt.subplots(figsize=size)
    mem = members('linear', test, short=(size == ps.HALF or test == 'ar'))   # fig. 6: L_z only, the long labels hid data
    cols = ps.ramp(len(mem))
    for (lab, key, val), col in zip(mem, cols):
        e = load(d, 'error', dt, key, val)
        if e is None:
            continue
        if not corrected and len(e) > 1:
            e = e[:-1]
        ax.semilogy(np.arange(len(e)), e, color=col, label=lab, **ps.STAR)
    ps.finish(ax, xlabel='KSP iteration', ylabel='relative error')
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.savefig(os.path.join(out, name)); plt.close(fig)


def linear_shift_sweep(d, corrected, out, name, size):
    fig, ax = plt.subplots(figsize=size)
    mem = members('linear', 'dt')
    cols = ps.ramp(len(mem))
    for (lab, key, c1), col in zip(mem, cols):
        xs, ys = [], []
        for dt in LIN['dts']:
            e = load(d, 'error', dt, 'shift', c1 / dt)
            if e is None:
                continue
            n = len(e) - 1
            its = n if corrected else len(e)
            if n >= MAXIT:
                continue
            xs.append(dt); ys.append(its)
        ax.plot(xs, ys, color=col, label=lab, **ps.STAR)
    ax.set_xscale('log')
    ps.finish(ax, xlabel=r'$\Delta t$ (s)', ylabel='KSP iterations', small=(size == ps.HALF))
    ax.legend(loc='upper right', ncol=3, columnspacing=1.0)   # fig. 7: three columns above the data
    ax.set_ylim(0, 120)
    fig.savefig(os.path.join(out, name)); plt.close(fig)


def nonlinear_history(d, test, dt, corrected, out, name, size, drop_last=False):
    '''drop_last: leave out the converged last Newton iterate (what plot_ar.py and the paper's AR dt = 5000 figure do).'''
    fig, ax = plt.subplots(figsize=size)
    mem = members('nonlinear', test, short=(size == ps.HALF))
    cols = ps.ramp(len(mem))
    for (lab, key, val), col in zip(mem, cols):
        tgt = val / dt if key == 'shift' else val
        se = load(d, 'snes_error', dt, key, tgt)
        cu = load(d, 'snes_ksp_cum', dt, key, tgt)
        if se is None or cu is None or se.size == 0:
            continue
        if (drop_last or corrected) and se.size > 1:
            se, cu = se[:-1], cu[:-1]
        ax.semilogy(cu, np.clip(se, 1e-16, None), color=col, label=lab, **ps.STAR)
    crowded = size == ps.HALF and len(mem) >= 6
    ps.finish(ax, xlabel='cumulative KSP iterations', ylabel=r'$\|U_k - U^*\| / \|U^*\|$', legend=not crowded, small=(size == ps.HALF))
    # crowded (fig. 11, seven C_1 in a half-width panel): no legend, the caption refers to the legend of fig. 12,
    # which uses the same colours for the same seven values (a legend in the panel hid data)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    fig.savefig(os.path.join(out, name)); plt.close(fig)


def nonlinear_shift_sweep(d, corrected, out, name, size):
    fig, ax = plt.subplots(figsize=size)
    mem = members('nonlinear', 'dt')
    cols = ps.ramp(len(mem))
    for (lab, key, c1), col in zip(mem, cols):
        xs, ys = [], []
        for dt in NL['dts']:
            e = load(d, 'error', dt, 'shift', c1 / dt)
            cu = load(d, 'snes_ksp_cum', dt, 'shift', c1 / dt)
            if e is None:
                continue
            if corrected:
                if cu is None:
                    continue
                per = np.diff(cu) if len(cu) > 1 else cu
                if np.any(per >= MAXIT):
                    continue
                its = int(cu[-1])
            else:
                its = len(e)
                if its >= MAXIT:
                    continue
            xs.append(dt); ys.append(its)
        ax.plot(xs, ys, color=col, label=lab, **ps.STAR)
    ax.set_xscale('log')
    ps.finish(ax, xlabel=r'$\Delta t$ (s)', ylabel='total KSP iterations' if corrected else 'KSP iterations', small=(size == ps.HALF))
    ax.set_ylim(bottom=0)
    fig.savefig(os.path.join(out, name)); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--set', choices=list(SUFFIX), default='paper')
    ap.add_argument('--out', default=None, help='output directory (default: ../../Patch-Preconditioner-for-Boussinesq/Image_restyled[/<set>])')
    ap.add_argument('--corrected', action='store_true', help='iteration counts n / total KSP and full histories; output to <out>/corrected')
    a = ap.parse_args()
    out = a.out or os.path.join(YR2, 'Patch-Preconditioner-for-Boussinesq', 'Image_restyled', '' if a.set == 'paper' else {'star': 'star_rerun', 'vanka': 'single_column'}[a.set])
    if a.corrected:
        out = os.path.join(out, 'corrected')
    os.makedirs(out, exist_ok=True)
    ps.apply()
    written = []

    def do(fn, *args):
        d = args[0]
        if not os.path.isdir(d):
            print('missing folder', d); return
        fn(*args)
        written.append([x for x in args if isinstance(x, str) and x.endswith('.png')][0])

    # linear (figs 4-7), placed at 0.8-0.85 linewidth
    do(linear_history, folder(a.set, 'linear', 'dx'), 'dx', 50.0, a.corrected, out, 'error_dx_Robust_dt50.0.png', ps.SINGLE)
    do(linear_history, folder(a.set, 'linear', 'dz'), 'dz', 50.0, a.corrected, out, 'error_dz_Robust_dt50.0.png', ps.SINGLE)
    do(linear_history, folder(a.set, 'linear', 'ar'), 'ar', 50.0, a.corrected, out, 'error_AR_Robust_dt50.0.png', ps.SINGLE)
    do(linear_shift_sweep, folder(a.set, 'linear', 'dt'), a.corrected, out, 'error_shift.png', ps.SINGLE)
    # nonlinear (figs 8-12), each panel at 0.48 linewidth
    do(nonlinear_history, folder(a.set, 'nonlinear', 'dx'), 'dx', 50.0, a.corrected, out, 'snes_error_dx_dt50.0.png', ps.HALF)
    do(nonlinear_history, folder(a.set, 'nonlinear', 'dz'), 'dz', 50.0, a.corrected, out, 'snes_error_dz_dt50.0.png', ps.HALF)
    do(nonlinear_history, folder(a.set, 'nonlinear', 'ar'), 'ar', 100.0, a.corrected, out, 'snes_error_AR_dt100.0.png', ps.HALF)
    do(nonlinear_history, folder(a.set, 'nonlinear', 'ar'), 'ar', 5000.0, a.corrected, out, 'snes_error_AR_dt5000.0.png', ps.HALF, True)   # the paper's figure drops the converged point
    do(nonlinear_history, folder(a.set, 'nonlinear', 'dt'), 'dt', 500.0, a.corrected, out, 'snes_error_shift_dt500.0.png', ps.HALF)
    do(nonlinear_shift_sweep, folder(a.set, 'nonlinear', 'dt'), a.corrected, out, 'error_shift_nonlinear.png', ps.HALF)
    print(f'{len(written)} figures written to {out}')


if __name__ == '__main__':
    main()
