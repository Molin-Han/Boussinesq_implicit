'''Re-draw each of the ten paper figures in the ORIGINAL (pre-restyle) matplotlib style, but from the data exactly as
restyle_paper_figures.py loads them (same folders, same file matching, same conventions), and diff against the
original PNG in Image_new/ or Image_nonlinear/. If the curves coincide pixel for pixel, the restyled figure shows
the same data as the original. Outputs go to this folder only.'''
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch')
import restyle_paper_figures as R

OUT = os.path.dirname(os.path.abspath(__file__)) + '/replica'
os.makedirs(OUT, exist_ok=True)
PAPER = '/home/mh2520/molin/phd/yr2/Patch-Preconditioner-for-Boussinesq'
L = 3.0e5

def _sci(x):
    e = int(np.floor(np.log10(abs(x)))); m = x / 10.0 ** e
    return rf'${m:.1f}\times 10^{{{e}}}$'

def lin_hist(test, labf, legfs, tight, name):
    plt.rcdefaults()
    d = R.folder('paper', 'linear', test)
    fig, ax = plt.subplots(figsize=(10, 6))
    for lab, key, val in R.members('linear', test):
        e = R.load(d, 'error', 50.0, key, val)
        e = e[:-1]
        ax.semilogy(np.arange(len(e)), e, marker='o', label=labf(val))
    ax.legend(fontsize=legfs); ax.set_ylabel('relative error', fontsize=20); ax.set_title('Error convergence at dt=50.0')
    fig.savefig(f'{OUT}/{name}', **({'bbox_inches': 'tight'} if tight else {})); plt.close(fig)

def shift_sweep(mode, name):
    plt.rcdefaults()
    S = R.LIN if mode == 'linear' else R.NL
    d = R.folder('paper', mode, 'dt')
    fig, ax = plt.subplots(figsize=(10, 6))
    for lab, key, c1 in R.members(mode, 'dt'):
        its = []
        for dt in S['dts']:
            e = R.load(d, 'error', dt, 'shift', c1 / dt)
            n = len(e)
            its.append(np.nan if n >= 100 else n)     # plot_shift.py: len(error), NaN if >= maxit (=100)
        ax.semilogx(S['dts'], its, marker='o', label=f'K={_sci(c1)}')
    ax.legend(fontsize=12); ax.set_ylabel('KSP iterations', fontsize=20); ax.set_title('KSP iterations vs dt for varying K')
    fig.savefig(f'{OUT}/{name}', bbox_inches='tight'); plt.close(fig)

def nl_old(test, dt, labf, name):
    '''older plot script (default figsize, x label, default legend, converged point kept)'''
    plt.rcdefaults()
    d = R.folder('paper', 'nonlinear', test)
    fig, ax = plt.subplots()
    for lab, key, val in R.members('nonlinear', test, short=True):
        tgt = val / dt if key == 'shift' else val
        se = R.load(d, 'snes_error', dt, key, tgt); cu = R.load(d, 'snes_ksp_cum', dt, key, tgt)
        ax.semilogy(cu, np.clip(se, 1e-16, None), marker='o', label=labf(val))
    ax.set_xlabel('cumulative KSP iterations'); ax.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$'); ax.legend()
    fig.savefig(f'{OUT}/{name}'); plt.close(fig)

def nl_ar5000(name):
    plt.rcdefaults()
    d = R.folder('paper', 'nonlinear', 'ar')
    fig, ax = plt.subplots(figsize=(10, 6))
    for lab, key, val in R.members('nonlinear', 'ar', short=True):
        se = R.load(d, 'snes_error', 5000.0, key, val)[:-1]; cu = R.load(d, 'snes_ksp_cum', 5000.0, key, val)[:-1]
        ax.semilogy(cu, np.clip(se, 1e-16, None), marker='o', label=f'AR={_sci(val)}')
    ax.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$', fontsize=20); ax.set_title('SNES error vs cumulative KSP iterations at dt=5000.0')
    ax.legend(fontsize=16); fig.savefig(f'{OUT}/{name}', bbox_inches='tight'); plt.close(fig)

jobs = [
    ('Image_new', 'error_dx_Robust_dt50.0.png', lambda n: lin_hist('dx', lambda v: f'nx={v}', 20, False, n)),
    ('Image_new', 'error_dz_Robust_dt50.0.png', lambda n: lin_hist('dz', lambda v: f'nz={v}', 20, False, n)),
    ('Image_new', 'error_AR_Robust_dt50.0.png', lambda n: lin_hist('ar', lambda v: f'AR={_sci(v)}', 16, True, n)),
    ('Image_new', 'error_shift.png', lambda n: shift_sweep('linear', n)),
    ('Image_nonlinear', 'snes_error_dx_dt50.0.png', lambda n: nl_old('dx', 50.0, lambda v: f'dx={L / v}', n)),
    ('Image_nonlinear', 'snes_error_dz_dt50.0.png', lambda n: nl_old('dz', 50.0, lambda v: f'dz={4000 / v}', n)),
    ('Image_nonlinear', 'snes_error_AR_dt100.0.png', lambda n: nl_old('ar', 100.0, lambda v: f'AR={round(v, 5)}', n)),
    ('Image_nonlinear', 'snes_error_AR_dt5000.0.png', lambda n: nl_ar5000(n)),
    ('Image_nonlinear', 'snes_error_shift_dt500.0.png', lambda n: nl_old('dt', 500.0, lambda v: f'C1={v}', n)),
    ('Image_nonlinear', 'error_shift_nonlinear.png', lambda n: shift_sweep('nonlinear', n)),
]
for folder, name, fn in jobs:
    fn(name)
    a = np.asarray(Image.open(f'{PAPER}/{folder}/{name}').convert('RGB')).astype(int)
    b = np.asarray(Image.open(f'{OUT}/{name}').convert('RGB')).astype(int)
    if a.shape != b.shape:
        print(f'{name:32s} SIZE DIFFERS original {a.shape} replica {b.shape}'); continue
    diff = np.abs(a - b).max(axis=2)
    nz = diff > 40
    Image.fromarray((255 - np.clip(diff * 3, 0, 255)).astype(np.uint8)).save(f'{OUT}/diff_{name}')
    ys, xs = np.nonzero(nz)
    box = (xs.min(), ys.min(), xs.max(), ys.max()) if nz.any() else None
    print(f'{name:32s} size {a.shape[1]}x{a.shape[0]}  max|diff|={diff.max():3d}  pixels with |diff|>40: {nz.sum():6d}  bbox={box}')
