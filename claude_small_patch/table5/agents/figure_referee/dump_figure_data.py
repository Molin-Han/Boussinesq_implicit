'''Print exactly the data that restyle_paper_figures.py (--set paper, default conventions) puts in each of the ten
figures, plus the true counts (n = lines-1 for linear; total KSP and Newton steps for nonlinear). Read-only.'''
import os, sys
import numpy as np
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch')
import restyle_paper_figures as R

def lin_hist(test):
    d = R.folder('paper', 'linear', test)
    print(f'== linear {test} dt=50  ({d})')
    for lab, key, val in R.members('linear', test):
        f = R.find(d, 'error', 50.0, key, val)
        e = np.atleast_1d(np.loadtxt(f))
        shown = e[:-1]
        print(f'  {lab:45s} file={os.path.basename(f):45s} lines={len(e):3d} n={len(e)-1:3d} curve ends at x={len(shown)-1:3d} '
              f'last shown={shown[-1]:.2e} dropped final={e[-1]:.1e}')

def lin_shift():
    d = R.folder('paper', 'linear', 'dt')
    print(f'== linear shift sweep ({d}); plotted = len(e) = n+1, omitted if n >= 100')
    for lab, key, c1 in R.members('linear', 'dt'):
        row = []
        for dt in R.LIN['dts']:
            f = R.find(d, 'error', dt, 'shift', c1 / dt)
            if f is None:
                row.append(f'{dt:g}:nofile'); continue
            e = np.atleast_1d(np.loadtxt(f)); n = len(e) - 1
            row.append(f'{dt:g}:{len(e) if n < 100 else "CAP("+str(len(e))+")"}')
        print(f'  {lab:20s}', ' '.join(row))

def nl_hist(test, dt, drop_last=False):
    d = R.folder('paper', 'nonlinear', test)
    print(f'== nonlinear {test} dt={dt} drop_last={drop_last} ({d})')
    for lab, key, val in R.members('nonlinear', test, short=True):
        tgt = val / dt if key == 'shift' else val
        se = R.load(d, 'snes_error', dt, key, tgt); cu = R.load(d, 'snes_ksp_cum', dt, key, tgt)
        if se is None or cu is None:
            print(f'  {lab}: MISSING'); continue
        fse = R.find(d, 'snes_error', dt, key, tgt)
        s_se, s_cu = (se[:-1], cu[:-1]) if drop_last and se.size > 1 else (se, cu)
        print(f'  {lab:28s} {os.path.basename(fse):40s} cum={list(map(int, cu))} err={[f"{x:.1e}" for x in se]} '
              f'shown_last=({int(s_cu[-1])},{s_se[-1]:.1e}) newton_steps={len(cu)-1}')

def nl_shift():
    d = R.folder('paper', 'nonlinear', 'dt')
    print(f'== nonlinear shift sweep ({d}); plotted = len(error) (= total KSP + Newton steps), omitted if >= 100')
    for lab, key, c1 in R.members('nonlinear', 'dt'):
        row = []
        for dt in R.NL['dts']:
            fe = R.find(d, 'error', dt, 'shift', c1 / dt)
            cu = R.load(d, 'snes_ksp_cum', dt, 'shift', c1 / dt)
            if fe is None:
                row.append(f'{dt:g}:nofile'); continue
            e = np.atleast_1d(np.loadtxt(fe))
            tot = None if cu is None else int(cu[-1]); ns = None if cu is None else len(cu) - 1
            per = None if cu is None else np.diff(cu).astype(int).tolist()
            flag = '' if len(e) < 100 else '(OMITTED)'
            row.append(f'{dt:g}:{len(e)}{flag}[tot={tot},N={ns},per={per}]')
        print(f'  {lab:22s}', ' '.join(row))

lin_hist('dx'); lin_hist('dz'); lin_hist('ar'); lin_shift()
nl_hist('dx', 50.0); nl_hist('dz', 50.0); nl_hist('ar', 100.0); nl_hist('ar', 5000.0, True); nl_hist('dt', 500.0); nl_shift()
