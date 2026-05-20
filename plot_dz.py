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

parser = ArgumentParser(
    description='Shifted simplified steady Linear Boussinesq equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--maxit', type=int, default=150, help='Max iteration number for the first ksp of the linear solve.')
parser.add_argument('--dts', nargs='+', type=float, default=1.0, help='The time step looping for.')
parser.add_argument('--nzs', nargs='+', type=int, default=1.0, help='The number of elements vertical-wise looping for.')
parser.add_argument('--C1', type=float, default=0.001, help='The default shift constant.')
parser.add_argument('--height', type=float, default=4000, help='The default height of the domain.')

args = parser.parse_known_args()
args = args[0]

dts = args.dts
nzs = args.nzs
C1 = args.C1
height = args.height

T = 1e4
dts_scaled = (np.array(dts)/T).tolist()

fig, ax = plt.subplots()
fig_scale, ax_scale = plt.subplots()

fig_res, ax_res = plt.subplots()
fig_res_scale, ax_res_scale = plt.subplots()

for nz in nzs:
    it_list = []
    it_res_list = []
    for dt in dts:
        deltaz = height / nz
        error = _load_or_zeros('error', dt, 'nz', nz)
        residual = _load_or_zeros('residual', dt, 'nz', nz)
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
    ax.semilogx(dts, it_list, marker='o', label=f'nz={nz}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
    ax_scale.semilogx(dts_scaled, it_list, marker='o', label=f'nz={nz}')
    ax_scale.legend()
    ax_scale.set_xlabel('dt')
    ax_scale.set_ylabel('its')

    ax_res.semilogx(dts, it_res_list, marker='o', label=f'nz={nz}')
    ax_res.legend()
    ax_res.set_xlabel('dt')
    ax_res.set_ylabel('its')
    ax_res_scale.semilogx(dts_scaled, it_res_list, marker='o', label=f'nz={nz}')
    ax_res_scale.legend()
    ax_res_scale.set_xlabel('dt')
    ax_res_scale.set_ylabel('its')
fig.savefig("error_dz.png")
fig_scale.savefig("error_dz_scaled_t.png")
fig_res.savefig("residual_dz.png")
fig_res_scale.savefig("residual_dz_scaled_t.png")

for dt in dts:
    it_list = []
    it_res_list = []
    fig_rob, ax_rob = plt.subplots()
    fig_res_rob, ax_res_rob = plt.subplots()
    for nz in nzs:
        deltaz = height / nz
        error = _load_or_zeros('error', dt, 'nz', nz)
        residual = _load_or_zeros('residual', dt, 'nz', nz)
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
        x = np.arange(its)
        x_res = np.arange(its_res)
        ax_rob.semilogy(x, error, label=f'nz={nz}')
        ax_rob.legend()
        ax_res_rob.semilogy(x_res, residual, label=f'nz={nz}')
        ax_res_rob.legend()
        plt.xlabel('its_num')
        plt.ylabel('log_error')
    fig_rob.savefig(f'error_dz_Robust_dt{dt}.png')
    fig_res_rob.savefig(f'residual_dz_Robust_dt{dt}.png')
    plt.close()


# SNES error vs cumulative KSP iteration count, one curve per nz value, per dt.
for dt in dts:
    fig_snes, ax_snes = plt.subplots()
    has_data = False
    for nz in nzs:
        deltaz = height / nz
        snes_err_path = _find_data_path('snes_error', dt, 'nz', nz)
        snes_ksp_path = _find_data_path('snes_ksp_cum', dt, 'nz', nz)
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
        snes_err = np.clip(snes_err, 1e-16, None)
        ax_snes.semilogy(snes_ksp_cum, snes_err, marker='o', label=f'dz={deltaz}')
        has_data = True
    if has_data:
        ax_snes.set_xlabel('cumulative KSP iterations')
        ax_snes.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$')
        ax_snes.legend()
        fig_snes.savefig(f'snes_error_dz_dt{dt}.png')
    plt.close(fig_snes)







