'''
Fig. 14(a) of the paper (aspect-ratio sweep, K = 1, total FGMRES iterations against dt, star and single-column
patch) with the committed preconditioner and with the corrected one (--pc_fix of Nonlinear_pcfix.py).
Reads ../../../data_claude_paper_nonlinear_ar_patch{star,vanka}{,_pcfix}/summary_iterations.csv, writes
fig14a_pcfix.png and fig14a_pcfix.md next to this file.
'''
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
YR2 = os.path.abspath(os.path.join(HERE, '..', '..', '..'))


def table(folder):
    path = f'{YR2}/{folder}/summary_iterations.csv'
    if not os.path.exists(path):
        return None
    out = {}
    for r in csv.DictReader(open(path)):
        tot = r['total_ksp']
        out[(float(r['dt']), float(r['height']))] = (int(float(tot)) if tot not in ('', 'nan') else None, r['status'], r['newton_steps'])
    return out


T = {(patch, fix): table(f'data_claude_paper_nonlinear_ar_patch{patch}{"_pcfix" if fix else ""}')
     for patch in ('star', 'vanka') for fix in (False, True)}
ref = T[('star', False)]
dts = sorted({k[0] for k in ref}); heights = sorted({k[1] for k in ref}, reverse=True)
RAMP = plt.get_cmap('viridis')
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5), sharey=True)
for ax, fix in zip(axs, (False, True)):
    for i, h in enumerate(heights):
        col = RAMP(i / (len(heights) - 1))
        for patch, ls, mfc in (('star', '-', col), ('vanka', '--', 'white')):
            tab = T[(patch, fix)]
            if tab is None:
                continue
            xs = [dt for dt in dts if tab.get((dt, h)) and tab[(dt, h)][1] == 'ok' and tab[(dt, h)][0]]
            ax.plot(xs, [tab[(dt, h)][0] for dt in xs], ls, marker='o', color=col, mfc=mfc, label=f'{h/1000:g} km' if patch == 'star' else None)
    ax.set_xscale('log'); ax.set_xlabel(r'$\Delta t$ (s)'); ax.grid(alpha=0.3)
    ax.set_title(('committed preconditioner (fig. 14(a))' if not fix else 'with the pressure back-substitution (--pc_fix)'), fontsize=10)
    ax.plot([], [], '-', marker='o', color='grey', label='star'); ax.plot([], [], '--', marker='o', mfc='white', color='grey', label='single column')
    ax.legend(fontsize=7, ncol=2)
axs[0].set_ylabel('total FGMRES iterations over the Newton steps')
fig.suptitle(r'Aspect-ratio sweep, $\delta = \delta_{thr}(\Delta t, \varepsilon, N)$ ($K = 1$), nonlinear code, no EW', fontsize=11)
fig.savefig(os.path.join(HERE, 'fig14a_pcfix.png'), dpi=160, bbox_inches='tight')

lines = ['# Fig. 14(a) data with and without the pressure back-substitution (total KSP / Newton steps)\n',
         'cell = star | column | star+fix | column+fix\n', '| dt | ' + ' | '.join(f'Lz={h/1000:g} km' for h in heights) + ' |', '|---|' + '---|' * len(heights)]
for dt in dts:
    cells = []
    for h in heights:
        parts = []
        for key in (('star', False), ('vanka', False), ('star', True), ('vanka', True)):
            tab = T[key]
            v = tab.get((dt, h)) if tab else None
            parts.append('-' if v is None else (str(v[0]) if v[1] == 'ok' and v[0] else v[1]))
        cells.append(' / '.join(parts))
    lines.append(f'| {dt:g} | ' + ' | '.join(cells) + ' |')
# ratios
for key, name in ((('star', True), 'star: fixed / committed'), (('vanka', True), 'column: fixed / committed'), (('vanka', False), 'column / star, committed'), (('vanka', True), 'column / star, fixed')):
    pass
def ratios(a, b):
    r = [a[k][0] / b[k][0] for k in b if a and k in a and a[k][1] == 'ok' and b[k][1] == 'ok' and a[k][0] and b[k][0]]
    return (f'n={len(r)}, median {np.median(r):.2f}, min {np.min(r):.2f}, max {np.max(r):.2f}') if r else 'no data'
lines.append('')
lines.append(f"star fixed/committed: {ratios(T[('star', True)], T[('star', False)])}")
lines.append(f"column fixed/committed: {ratios(T[('vanka', True)], T[('vanka', False)])}")
lines.append(f"column/star committed: {ratios(T[('vanka', False)], T[('star', False)])}")
lines.append(f"column/star fixed: {ratios(T[('vanka', True)], T[('star', True)])}")
open(os.path.join(HERE, 'fig14a_pcfix.md'), 'w').write('\n'.join(lines) + '\n')
print('written fig14a_pcfix.png / .md')
