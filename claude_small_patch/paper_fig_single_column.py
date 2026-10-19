'''
Figure for the paper's single-column section: the paper's own robustness sweeps (delta = C_1/dt) for the
star patch (solid, filled) and the single-column patch (dashed, open), both codes, three panels each:
aspect ratio, mesh resolution at dt = 50 s, shift sweep.  Writes
../../Patch-Preconditioner-for-Boussinesq/Image_restyled/single_column_{linear,nonlinear}.png (1 x 3 panels each).

Data: ../../data_<mode>_<test>_deg2[_ew]_{star,vanka} (the 2026-10-16 re-runs with identical flags), read
through compare_paper_patches.record() so that the conventions (counts, capped runs) are the same as in
claude_small_patch/results/.
'''
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import paper_style as ps                                   # noqa: E402
import compare_paper_patches as cpp                        # noqa: E402

OUT = os.path.join(HERE, '..', '..', 'Patch-Preconditioner-for-Boussinesq', 'Image_restyled', 'single_column_sweeps.png')
ps.apply()


def table(mode, test):
    S = cpp.SWEEP[mode]
    mem = cpp.members(mode, test)
    data = {}
    for patch in ('star', 'vanka'):
        d = cpp.folder(mode, test, patch)
        fj = cpp.failed_jobs(d)
        data[patch] = {(val, dt): cpp.record(mode, test, patch, key, val, dt, fj) for (_, key, val) in mem for dt in S['dts']}
    return mem, data


LEG = dict(fontsize=6.5, handlelength=1.6, borderpad=0.3, labelspacing=0.25, columnspacing=0.8)
OUTDIR = os.path.dirname(OUT)
for mode in ('linear', 'nonlinear'):
    S = cpp.SWEEP[mode]
    ylab = 'KSP iterations' if mode == 'linear' else 'total KSP iterations'
    fig, axes = plt.subplots(1, 3, figsize=(6.8, 2.6))
    # (a) aspect ratio: label the thickest, a middle and the thinnest domain only
    ax = axes[0]
    mem, data = table(mode, 'ar')
    cols = ps.ramp(len(mem))
    show = {0, len(mem) // 2, len(mem) - 1}
    for i, (lab, _, val) in enumerate(mem):
        for patch, style in (('star', ps.STAR), ('vanka', ps.COLUMN)):
            xs = [dt for dt in S['dts'] if data[patch][(val, dt)]['status'] == 'ok']
            ax.plot(xs, [data[patch][(val, dt)]['count'] for dt in xs], color=cols[i],
                    label=lab if (patch == 'star' and i in show) else None, **style)
    ps.patch_legend(ax, star_label='star', column_label='single column')
    ax.set_xscale('log'); ax.set_title(f"(a) aspect ratio, $C_1 = {ps.sci(S['C1'])[1:-1]}$")
    ax.set_ylim(0, 1.75 * max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(ncol=1, loc='upper right', **LEG)
    ax.set_xlabel(r'$\Delta t$ (s)'); ax.set_ylabel(ylab)
    # (b) mesh at dt = 50
    ax = axes[1]
    for test, mk, colr, name in (('dx', 'o', 'C0', '$N_x$'), ('dz', 's', 'C1', '$N_z$')):
        mem, data = table(mode, test)
        ns = [val for (_, _, val) in mem]
        for patch, ls, mfc, pname in (('star', '-', colr, 'star'), ('vanka', '--', 'white', 'single column')):
            ax.plot(ns, [data[patch][(n, 50.0)]['count'] for n in ns], ls=ls, marker=mk, color=colr, mfc=mfc, label=f'{pname}, {name}')
    ax.set_title(r"(b) mesh resolution, $\Delta t = 50$ s")
    ax.set_xticks([10, 20, 40, 60, 80]); ax.set_xlabel(r'$N_x$ (coarse mesh) or $N_z$')
    ax.set_ylim(0, 1.6 * max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(loc='upper left', ncol=2, **LEG)
    # (c) shift sweep: label the largest, a middle and the smallest C_1 only
    ax = axes[2]
    mem, data = table(mode, 'dt')
    cols = ps.ramp(len(mem))
    show = {0, len(mem) // 2, len(mem) - 1}
    for i, (lab, _, val) in enumerate(mem):
        for patch, style in (('star', ps.STAR), ('vanka', ps.COLUMN)):
            xs = [dt for dt in S['dts'] if data[patch][(val, dt)]['status'] == 'ok']
            ax.plot(xs, [data[patch][(val, dt)]['count'] for dt in xs], color=cols[i],
                    label=lab if (patch == 'star' and i in show) else None, **style)
    ax.set_xscale('log'); ax.set_title(r"(c) shift $\delta = C_1/\Delta t$")
    ax.set_ylim(0, 1.45 * max(max(l.get_ydata()) for l in ax.get_lines() if len(l.get_ydata())))
    ax.legend(ncol=1, loc='upper right', **LEG)
    ax.set_xlabel(r'$\Delta t$ (s)')
    for a_ in axes:
        a_.grid(True, alpha=0.3)
    fig.tight_layout(w_pad=0.8)
    out = os.path.join(OUTDIR, f'single_column_{mode}.png')
    fig.savefig(out); plt.close(fig)
    print('written', out)
