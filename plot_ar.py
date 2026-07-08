from firedrake import *
import glob
import re
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from firedrake.output import VTKFile
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


def _find_data_path(prefix, dt, key, target, ext='.out'):
    # Pick the file {prefix}_dt{dt}_{key}<value>{ext} whose <value> is closest
    # in log-space to target. Returns None if no candidate within ~20% (log).
    candidates = glob.glob(f'{prefix}_dt{dt}_{key}*{ext}')
    if not candidates:
        return None
    rx = re.compile(rf'_{re.escape(key)}([^/]+){re.escape(ext)}$')
    best, best_d = None, float('inf')
    for p in candidates:
        m = rx.search(p)
        if not m:
            continue
        try:
            v = float(m.group(1))
        except ValueError:
            continue
        if v > 0 and target > 0:
            d = abs(np.log(v) - np.log(target))
        else:
            d = abs(v - target)
        if d < best_d:
            best, best_d = p, d
    return best if best is not None and best_d <= 0.2 else None


def _load_or_zeros(prefix, dt, key, target, fallback_len=5):
    path = _find_data_path(prefix, dt, key, target)
    if path is None:
        return np.zeros(fallback_len)
    try:
        return np.loadtxt(path)
    except (OSError, ValueError):
        return np.zeros(fallback_len)


def _sci(x):
    if x == 0:
        return r'$0$'
    exp = int(np.floor(np.log10(abs(x))))
    mantissa = x / 10.0**exp
    return rf'${mantissa:.1f}\times 10^{{{exp}}}$'

parser = ArgumentParser(
    description='Shifted simplified steady Linear Boussinesq equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--maxit', type=int, default=100, help='Max iteration number for the first ksp of the linear solve.')
parser.add_argument('--dts', nargs='+', type=float, default=1.0, help='The time step looping for.')
parser.add_argument('--heights', nargs='+', type=float, default=1.0, help='The height looping for.')
parser.add_argument('--C1', type=float, default=0.001, help='The default shift constant.')
parser.add_argument('--length', type=float, default=3.0e5, help='The default length of the domain.')

args = parser.parse_known_args()
args = args[0]

dts = args.dts
heights = args.heights
C1 = args.C1
length = args.length

T = 1e4
dts_scaled = (np.array(dts)/T).tolist()

fig, ax = plt.subplots(figsize=(10, 6))
fig_scale, ax_scale = plt.subplots(figsize=(10, 6))

fig_res, ax_res = plt.subplots(figsize=(10, 6))
fig_res_scale, ax_res_scale = plt.subplots(figsize=(10, 6))

for height in heights:
    it_list = []
    it_res_list = []
    i = 0
    for dt in dts:
        i += 1
        ar = height / length
        error = _load_or_zeros('error', dt, 'ar', ar)
        residual = _load_or_zeros('residual', dt, 'ar', ar)
        its = len(error)
        its_res = len(residual)
        if its >= args.maxit:
            it_list.append(np.nan)
            # it_list.append(its)
        else:
            it_list.append(its)
        if its_res >= args.maxit:
            it_res_list.append(np.nan)
            # it_list.append(its)
        else:
            it_res_list.append(its_res)
    ax.semilogx(dts, it_list, marker='o', label=f'AR={_sci(ar)}')
    ax.legend(fontsize=16)
    # ax.set_xlabel('dt')
    ax.set_ylabel('KSP iterations', fontsize=20)
    ax.set_title('KSP iterations vs dt for varying aspect ratio')
    ax_scale.semilogx(dts_scaled, it_list, marker='o', label=f'AR={_sci(ar)}')
    ax_scale.legend(fontsize=16)
    # ax_scale.set_xlabel('dt / T')
    ax_scale.set_ylabel('KSP iterations', fontsize=20)
    ax_scale.set_title('KSP iterations vs dt/T for varying aspect ratio')
    ax_res.semilogx(dts, it_res_list, marker='o', label=f'AR={_sci(ar)}')
    ax_res.legend(fontsize=16)
    # ax_res.set_xlabel('dt')
    ax_res.set_ylabel('KSP iterations (residual)', fontsize=20)
    ax_res.set_title('Residual KSP iterations vs dt for varying aspect ratio')
    ax_res_scale.semilogx(dts_scaled, it_res_list, marker='o', label=f'AR={_sci(ar)}')
    ax_res_scale.legend(fontsize=16)
    # ax_res_scale.set_xlabel('dt / T')
    ax_res_scale.set_ylabel('KSP iterations (residual)', fontsize=20)
    ax_res_scale.set_title('Residual KSP iterations vs dt/T for varying aspect ratio')
fig.savefig("error_AR.png", bbox_inches='tight')
fig_scale.savefig("error_AR_scaled_t.png", bbox_inches='tight')
fig_res.savefig("res_AR.png", bbox_inches='tight')
fig_res_scale.savefig("res_AR_scaled_t.png", bbox_inches='tight')

for dt in dts:
    it_list = []
    it_res_list = []
    fig_rob, ax_rob = plt.subplots(figsize=(10, 6))
    fig_res_rob, ax_res_rob = plt.subplots(figsize=(10, 6))
    for height in heights:
        ar = height / length
        error = _load_or_zeros('error', dt, 'ar', ar)
        residual = _load_or_zeros('residual', dt, 'ar', ar)
        its = len(error)
        its_res = len(residual)
        if its >= args.maxit:
            it_list.append(np.nan)
        else:
            it_list.append(its)
        if its_res >= args.maxit:
            it_res_list.append(np.nan)
        else:
            it_res_list.append(its_res)
        err_plot = error[:-1] if len(error) > 1 else error
        res_plot = residual[:-1] if len(residual) > 1 else residual
        x = np.arange(len(err_plot))
        x_res = np.arange(len(res_plot))
        ax_rob.semilogy(x, err_plot, marker='o', label=f'AR={_sci(ar)}')
        ax_res_rob.semilogy(x_res, res_plot, marker='o', label=f'AR={_sci(ar)}')
    ax_rob.legend(fontsize=16)
    # ax_rob.set_xlabel('KSP iteration')
    ax_rob.set_ylabel('relative error', fontsize=20)
    ax_rob.set_title(f'Error convergence at dt={dt}')
    ax_res_rob.legend(fontsize=16)
    # ax_res_rob.set_xlabel('KSP iteration')
    ax_res_rob.set_ylabel('residual', fontsize=20)
    ax_res_rob.set_title(f'Residual convergence at dt={dt}')
    fig_rob.savefig(f'error_AR_Robust_dt{dt}.png', bbox_inches='tight')
    fig_res_rob.savefig(f'residual_AR_Robust_dt{dt}.png', bbox_inches='tight')
    plt.close(fig_rob)
    plt.close(fig_res_rob)


# SNES error vs cumulative KSP iteration count, one curve per AR value, per dt.
for dt in dts:
    fig_snes, ax_snes = plt.subplots(figsize=(10, 6))
    has_data = False
    for height in heights:
        ar = height / length
        snes_err_path = _find_data_path('snes_error', dt, 'ar', ar)
        snes_ksp_path = _find_data_path('snes_ksp_cum', dt, 'ar', ar)
        if snes_err_path is None or snes_ksp_path is None:
            continue
        try:
            snes_err = np.loadtxt(snes_err_path)
            snes_ksp_cum = np.loadtxt(snes_ksp_path)
        except (OSError, ValueError):
            continue
        snes_err = np.atleast_1d(snes_err)
        snes_ksp_cum = np.atleast_1d(snes_ksp_cum)
        if snes_err.size == 0:
            continue
        if snes_err.size > 1:
            snes_err = snes_err[:-1]
            snes_ksp_cum = snes_ksp_cum[:-1]
        snes_err = np.clip(snes_err, 1e-16, None)
        ax_snes.semilogy(snes_ksp_cum, snes_err, marker='o',
                         label=f'AR={_sci(ar)}')
        has_data = True
    if has_data:
        # ax_snes.set_xlabel('cumulative KSP iterations')
        ax_snes.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$', fontsize=20)
        ax_snes.set_title(f'SNES error vs cumulative KSP iterations at dt={dt}')
        ax_snes.legend(fontsize=16)
        fig_snes.savefig(f'snes_error_AR_dt{dt}.png', bbox_inches='tight')
    plt.close(fig_snes)


