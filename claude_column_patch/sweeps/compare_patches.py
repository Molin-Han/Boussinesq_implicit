'''
Star patch versus single-column (vanka) patch: tables and figures from the sweeps

    ../data_claude_paper_nonlinear_<test>_patchstar     (ASMStarPC,  construct_dim 0)
    ../data_claude_paper_nonlinear_<test>_patchvanka    (ASMVankaPC, construct_dim 1)

made by run_claude_paper.sh (-s _patchstar / -s _patchvanka).  Optionally the earlier
star sweeps ../data_claude_paper_nonlinear_<test> are read as a regression check of the
--patch edit (the star runs must reproduce them).

    python compare_patches.py [--data_root DIR] [--tests dx dz ar dt] [--out DIR]

Everything is read from summary_iterations.csv, which plot_claude_paper.py writes into
each data folder (total_ksp = total number of FGMRES iterations over all Newton steps of
the monitored time step; status = ok | capped | diverged | crash | missing).
'''
import os
import csv
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('--data_root', default=os.path.abspath(os.path.join(here, '..', '..', '..')), help='Folder holding the data_claude_paper_* folders.')
parser.add_argument('--tests', nargs='+', default=['dx', 'dz', 'ar', 'dt'])
parser.add_argument('--out', default=here, help='Where the tables and figures are written.')
args = parser.parse_args()

MEMBER_LABEL = {'dx': 'N_x (coarse)', 'dz': 'N_z', 'ar': 'L_z (m)', 'dt': 'K'}


def read(folder):
    path = os.path.join(folder, 'summary_iterations.csv')
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return list(csv.DictReader(f))


def member_value(test, row):
    return {'dx': int(row['nx']), 'dz': int(row['nz']), 'ar': float(row['height']), 'dt': float(row['K'])}[test]


def table(rows, test):
    '''{(dt, member): (total_ksp or None, newton_steps, status)}'''
    out = {}
    for r in rows:
        tot = r['total_ksp']
        out[(float(r['dt']), member_value(test, r))] = (
            int(float(tot)) if tot not in ('', 'nan') else None,
            int(float(r['newton_steps'])) if r['newton_steps'] not in ('', 'nan') else None,
            r['status'])
    return out


def cell(entry):
    if entry is None:
        return 'n/a'
    tot, _, status = entry
    if status == 'ok':
        return str(tot)
    return {'capped': '>100', 'diverged': 'div'}.get(status, status)


md = ['# Star patch vs single-column (vanka) patch: iteration counts\n',
      'Total FGMRES iterations over all Newton steps of the monitored time step, `star / vanka`.  '
      '`>100` = a linear solve hit the cap of 100 iterations, `div` = Newton diverged.  '
      'Shift $\\delta = K\\,\\delta_{thr}$; $K = 1$ in the dx, dz and ar sweeps.\n']
summary = []
for test in args.tests:
    star = read(os.path.join(args.data_root, f'data_claude_paper_nonlinear_{test}_patchstar'))
    vanka = read(os.path.join(args.data_root, f'data_claude_paper_nonlinear_{test}_patchvanka'))
    old = read(os.path.join(args.data_root, f'data_claude_paper_nonlinear_{test}'))
    if star is None or vanka is None:
        print(f'{test}: summary_iterations.csv missing for star or vanka, skipped')
        continue
    S, V = table(star, test), table(vanka, test)
    dts = sorted({k[0] for k in S} | {k[0] for k in V})
    members = sorted({k[1] for k in S} | {k[1] for k in V})
    if test == 'ar':
        members = members[::-1]

    md.append(f'\n## {test} sweep\n')
    md.append('| dt | ' + ' | '.join(f'{MEMBER_LABEL[test]} = {m:g}' for m in members) + ' |')
    md.append('|---|' + '---|' * len(members))
    ratios = []
    for dt in dts:
        cells = []
        for m in members:
            s, v = S.get((dt, m)), V.get((dt, m))
            cells.append(f'{cell(s)} / {cell(v)}')
            if s and v and s[2] == 'ok' and v[2] == 'ok' and s[0]:
                ratios.append((dt, m, s[0], v[0], v[0] / s[0]))
        md.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')
    if ratios:
        r = np.array([x[4] for x in ratios])
        worst = max(ratios, key=lambda x: x[4])
        n_s_fail = sum(1 for k in S if S[k][2] != 'ok')
        n_v_fail = sum(1 for k in V if V[k][2] != 'ok')
        md.append(f'\nRatio vanka/star over the {len(ratios)} runs where both converged: median {np.median(r):.2f}, '
                  f'mean {r.mean():.2f}, max {r.max():.2f} (dt = {worst[0]:g}, {MEMBER_LABEL[test]} = {worst[1]:g}: '
                  f'{worst[2]} -> {worst[3]}).  Not converged: star {n_s_fail}, vanka {n_v_fail} of {len(S)}.')
        summary.append((test, len(ratios), np.median(r), r.mean(), r.max(), n_s_fail, n_v_fail, len(S)))

    # ! Regression of the --patch edit: the new star sweep against the earlier star sweep.
    if old is not None:
        O = table(old, test)
        common = [k for k in S if k in O]
        diff = [(k, O[k], S[k]) for k in common if (O[k][0], O[k][2]) != (S[k][0], S[k][2])]
        md.append(f'\nRegression against the earlier star sweep (`data_claude_paper_nonlinear_{test}`): '
                  f'{len(common) - len(diff)} of {len(common)} common runs identical in total iterations and status.')
        for k, o, s in diff:
            md.append(f'  - dt = {k[0]:g}, {MEMBER_LABEL[test]} = {k[1]:g}: earlier {cell(o)}, now {cell(s)}')

    # ! Figure: total iterations against dt, one colour per member, solid = star, dashed = vanka.
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.8))
    ramp = ['#86b6ef', '#5598e7', '#2a78d6', '#256abf', '#1c5cab', '#104281', '#0d366b']   # one hue, light -> dark: ordered members
    for i, m in enumerate(members):
        col = ramp[round(i * (len(ramp) - 1) / max(1, len(members) - 1))]
        for P, ls, mk, name in ((S, '-', 'o', 'star'), (V, '--', 's', 'vanka')):
            xs = [dt for dt in dts if P.get((dt, m)) and P[(dt, m)][2] == 'ok']
            ys = [P[(dt, m)][0] for dt in xs]
            ax[0].plot(xs, ys, ls=ls, lw=2, marker=mk, ms=5, mfc=('white' if name == 'vanka' else col), color=col, label=f'{name}, {MEMBER_LABEL[test]} = {m:g}')
        xs = [dt for dt in dts if S.get((dt, m)) and V.get((dt, m)) and S[(dt, m)][2] == 'ok' and V[(dt, m)][2] == 'ok' and S[(dt, m)][0]]
        ax[1].plot(xs, [V[(dt, m)][0] / S[(dt, m)][0] for dt in xs], lw=2, marker='o', ms=5, color=col, label=f'{MEMBER_LABEL[test]} = {m:g}')
    ax[0].set_xscale('log'); ax[0].set_xlabel(r'$\Delta t$ (s)'); ax[0].set_ylabel('total FGMRES iterations')
    ax[0].set_title(f'{test} sweep: star (solid, filled) and vanka (dashed, open)'); ax[0].grid(alpha=0.3)
    ax[0].set_ylim(bottom=0)
    ax[0].legend(fontsize=7, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.16), frameon=False)   # below the axes: never on top of the data
    ax[1].axhline(1.0, color='#6b6a63', lw=1.0)
    ax[1].set_ylim(bottom=0)
    ax[1].set_xscale('log'); ax[1].set_xlabel(r'$\Delta t$ (s)'); ax[1].set_ylabel('iterations vanka / star')
    ax[1].set_title('ratio'); ax[1].grid(alpha=0.3)
    ax[1].legend(fontsize=7, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.16), frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, f'compare_{test}.png'), dpi=150)
    plt.close(fig)
    md.append(f'\n![{test}](compare_{test}.png)\n')

# ! Collapse test of the analysis: the penalty of the single-column patch against
# ! rho_1 = m L_z^2 / (pi^2 dx^2), the weight of the H1-in-x part of the energy of the divergence-free
# ! fields relative to its L2-in-x part, for the gravest vertical mode (dx = finest-level column width).
N_FREQ = 1.0e-2
pts = []
for test in ('dx', 'dz', 'ar'):
    star = read(os.path.join(args.data_root, f'data_claude_paper_nonlinear_{test}_patchstar'))
    vanka = read(os.path.join(args.data_root, f'data_claude_paper_nonlinear_{test}_patchvanka'))
    if star is None or vanka is None:
        continue
    S, V = table(star, test), table(vanka, test)
    for r in star:
        k = (float(r['dt']), member_value(test, r))
        if k in V and S[k][2] == 'ok' and V[k][2] == 'ok' and S[k][0]:
            dt, Lx, Lz, nx = float(r['dt']), float(r['length']), float(r['height']), int(r['nx'])
            m = 1.0 + (dt * N_FREQ) ** 2 / 4.0
            dxf = Lx / (4 * nx)
            pts.append((test, dt, m * Lz ** 2 / (np.pi ** 2 * dxf ** 2), S[k][0], V[k][0]))
if pts:
    colours = {'dx': '#2a78d6', 'dz': '#eb6834', 'ar': '#1baf7a'}
    markers = {'dx': 'o', 'dz': 's', 'ar': '^'}
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    for test in ('dx', 'dz', 'ar'):
        q = [p for p in pts if p[0] == test]
        ax[0].scatter([p[2] for p in q], [p[4] / p[3] for p in q], s=36, marker=markers[test], color=colours[test],
                      edgecolor='white', linewidth=0.8, label=f'{test} sweep ({len(q)} runs)')
        ax[1].scatter([p[2] for p in q], [p[4] - p[3] for p in q], s=36, marker=markers[test], color=colours[test],
                      edgecolor='white', linewidth=0.8, label=f'{test} sweep')
    for a in ax:
        a.set_xscale('log'); a.grid(alpha=0.3); a.legend(fontsize=8)
        a.set_xlabel(r'$\rho_1 = m L_z^2 / (\pi^2 \Delta x^2)$')
    ax[0].axhline(1.0, color='#6b6a63', lw=1.0); ax[0].set_ylim(bottom=0)
    ax[0].set_ylabel('iterations vanka / star'); ax[0].set_title('K = 1: ratio of total FGMRES iterations')
    ax[1].axhline(0.0, color='#6b6a63', lw=1.0)
    ax[1].set_ylabel('iterations vanka - star'); ax[1].set_title('K = 1: extra iterations of the single-column patch')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'collapse_rho1.png'), dpi=150)
    plt.close(fig)
    rho = np.array([p[2] for p in pts]); ratio = np.array([p[4] / p[3] for p in pts]); extra = np.array([p[4] - p[3] for p in pts])
    md.append('\n## Penalty of the single-column patch against $\\rho_1 = m L_z^2/(\\pi^2\\Delta x^2)$ (all $K = 1$ runs)\n')
    md.append('| $\\rho_1$ | runs | median vanka/star | max vanka/star | median extra iterations | max extra |')
    md.append('|---|---|---|---|---|---|')
    edges = [0, 1e-2, 1e-1, 1, 10, 100, np.inf]
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (rho >= lo) & (rho < hi)
        if sel.any():
            md.append(f'| [{lo:g}, {hi:g}) | {sel.sum()} | {np.median(ratio[sel]):.2f} | {ratio[sel].max():.2f} | {np.median(extra[sel]):.0f} | {extra[sel].max():.0f} |')
    md.append('\n![collapse](collapse_rho1.png)\n')

if summary:
    md.insert(2, '\n| sweep | runs compared | median vanka/star | mean | max | not converged (star) | not converged (vanka) |\n|---|---|---|---|---|---|---|\n'
              + '\n'.join(f'| {t} | {n} | {med:.2f} | {mean:.2f} | {mx:.2f} | {sf} of {tot} | {vf} of {tot} |'
                          for t, n, med, mean, mx, sf, vf, tot in summary) + '\n')
open(os.path.join(args.out, 'compare_iterations.md'), 'w').write('\n'.join(md) + '\n')
print('\n'.join(md))
