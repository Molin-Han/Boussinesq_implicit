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
parser.add_argument('--nzs', nargs='+', type=float, default=1.0, help='The number of elements vertical-wise looping for.')
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

for nz in nzs:
    it_list = []
    for dt in dts:
        deltaz = height / nz
        error = np.loadtxt(f'error_dt{dt}_nz{nz}.out')
        its = len(error)
        if its >= args.maxit:
            it_list.append(np.nan)
            # it_list.append(its)
        else:
            it_list.append(its)
    ax.semilogx(dts, it_list, label=f'dz={deltaz}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
    ax_scale.semilogx(dts_scaled, it_list, label=f'dz={deltaz}')
    ax_scale.legend()
    ax_scale.set_xlabel('dt')
    ax_scale.set_ylabel('its')
fig.savefig("error_dz.png")
fig_scale.savefig("error_dz_scaled_t.png")


for dt in dts:
    it_list = []
    fig_rob, ax_rob = plt.subplots()
    for nz in nzs:
        deltaz = height / nz
        error = np.loadtxt(f'error_dt{dt}_nz{nz}.out')
        its = len(error)
        if its >= args.maxit:
            it_list.append(np.nan)
            # it_list.append(its)
        else:
            it_list.append(its)
        x = np.arange(its)
        ax_rob.semilogy(x, error, label=f'nz={nz}')
        ax_rob.legend()
        plt.xlabel('its_num')
        plt.ylabel('log_error')
    fig_rob.savefig(f'error_dz_Robust_dt{dt}.png')
    plt.close()







