from firedrake import *
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from firedrake.output import VTKFile
from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter

parser = ArgumentParser(
    description='Shifted simplified steady Linear Boussinesq equation.',
    formatter_class=ArgumentDefaultsHelpFormatter
)

parser.add_argument('--maxit', type=int, default=150, help='Max iteration number for the first ksp of the linear solve.')
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

fig, ax = plt.subplots()
fig_scale, ax_scale = plt.subplots()

fig_res, ax_res = plt.subplots()
fig_res_scale, ax_res_scale = plt.subplots()

for height in heights:
    it_list = []
    it_res_list = []
    i = 0
    for dt in dts:
        i += 1
        ar = height / length
        try:
            error = np.loadtxt(f'error_dt{dt}_ar{ar}.out')
            residual = np.loadtxt(f'residual_dt{dt}_ar{ar}.out')
        except FileNotFoundError:
            error = np.zeros(5)
            residual = np.zeros(5)
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
    ax.semilogx(dts, it_list, label=f'AR={np.round(ar, decimals=5)}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
    ax_scale.semilogx(dts_scaled, it_list, label=f'AR={np.round(ar, decimals=5)}')
    ax_scale.legend()
    ax_scale.set_xlabel('dt')
    ax_scale.set_ylabel('its')
    ax_res.semilogx(dts, it_list, label=f'AR={np.round(ar, decimals=5)}')
    ax_res.legend()
    ax_res.set_xlabel('dt')
    ax_res.set_ylabel('its')
    ax_res_scale.semilogx(dts_scaled, it_list, label=f'AR={np.round(ar, decimals=5)}')
    ax_res_scale.legend()
    ax_res_scale.set_xlabel('dt')
    ax_res_scale.set_ylabel('its')
fig.savefig("error_AR.png")
fig_scale.savefig("error_AR_scaled_t.png")
fig_res.savefig("res_AR.png")
fig_res_scale.savefig("res_AR_scaled_t.png")

for dt in dts:
    it_list = []
    res_list = []
    fig_rob, ax_rob = plt.subplots()
    fig_res_rob, ax_res_rob = plt.subplots()
    for height in heights:
        ar = height / length
        try:
            error = np.loadtxt(f'error_dt{dt}_ar{ar}.out')
            residual = np.loadtxt(f'residual_dt{dt}_ar{ar}.out')
        except FileNotFoundError:
            error = np.zeros(5)
            residual = np.zeros(5)
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
        ax_rob.semilogy(x, error, label=f'AR={np.round(ar, decimals=5)}')
        ax_rob.legend()
        ax_res_rob.semilogy(x, residual, label=f'AR={np.round(ar, decimals=5)}')
        ax_res_rob.legend()
        plt.xlabel('its_num')
        plt.ylabel('log_error')
    fig_rob.savefig(f'error_AR_Robust_dt{dt}.png')
    fig_res_rob.savefig(f'residual_AR_Robust_dt{dt}.png')
    plt.close()


# SNES error vs cumulative KSP iteration count, one curve per AR value, per dt.
for dt in dts:
    fig_snes, ax_snes = plt.subplots()
    has_data = False
    for height in heights:
        ar = height / length
        try:
            snes_err = np.loadtxt(f'snes_error_dt{dt}_ar{ar}.out')
            snes_ksp_cum = np.loadtxt(f'snes_ksp_cum_dt{dt}_ar{ar}.out')
        except FileNotFoundError:
            continue
        snes_err = np.atleast_1d(snes_err)
        snes_ksp_cum = np.atleast_1d(snes_ksp_cum)
        if snes_err.size == 0:
            continue
        snes_err = np.clip(snes_err, 1e-16, None)
        ax_snes.semilogy(snes_ksp_cum, snes_err, marker='o',
                         label=f'AR={np.round(ar, decimals=5)}')
        has_data = True
    if has_data:
        ax_snes.set_xlabel('cumulative KSP iterations')
        ax_snes.set_ylabel(r'$\|U_k - U^*\| / \|U^*\|$')
        ax_snes.legend()
        fig_snes.savefig(f'snes_error_AR_dt{dt}.png')
    plt.close(fig_snes)


