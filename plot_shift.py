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
parser.add_argument('--C1_list',nargs='+', type=float, default=1.0, help='The C1 shift constant looping for.')

args = parser.parse_known_args()
args = args[0]

dts = args.dts
C1_list = args.C1_list

T = 1e4
print('!!!!',dts)
dts_scaled = (np.array(dts)/T).tolist()

fig, ax = plt.subplots()
fig_scale, ax_scale = plt.subplots()

for C1 in C1_list:
    it_list = []
    fig_rob, ax_rob = plt.subplots()
    for dt in dts:
        shift = np.round(C1 * dt**(-1.5), decimals=16)
        error = np.loadtxt(f'error_dt{dt}_shift{shift}.out')
        its = len(error)
        if its >= args.maxit:
            it_list.append(np.nan)
            # it_list.append(its)
        else:
            it_list.append(its)
        x = np.arange(its)
        ax_rob.semilogy(x, error, label=f'dt={dt}')
        ax_rob.legend()
        plt.xlabel('its_num')
        plt.ylabel('log_error')
    fig_rob.savefig(f'error_Robust_C1_{C1}.png')
    ax.semilogx(dts, it_list, label=f'C1={C1}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
    ax_scale.semilogx(dts_scaled, it_list, label=f'C1={C1}')
    ax_scale.legend()
    ax_scale.set_xlabel('dt')
    ax_scale.set_ylabel('its')
fig.savefig("error_shift.png")
fig_scale.savefig("error_shift_scaled_t.png")








