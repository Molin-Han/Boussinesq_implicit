from firedrake import *
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from firedrake.output import VTKFile

dts = [1.0, 5.0, 10.0, 50.0, 100.0, 200.0, 1000.0]
# C1_list = [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001]
# C1_list = [0.01, 0.001, 0.0001, 0.00001]
C1_list = [0.001]
fig, ax = plt.subplots()

for C1 in C1_list:
    it_list = []
    fig_tem, ax_tem = plt.subplots()
    fig_rob, ax_rob = plt.subplots()
    for dt in dts:
        if (C1 == 0.01 or C1== 0.1 or C1==1.0 or C1 == 0.001) and dt == 1.0:
            it_list.append(np.nan)
        else:
            shift = np.round(C1 * dt**(-1.5), decimals=10)
            error = np.loadtxt(f'error_dt{dt}_shift{shift}.out')
            its = len(error)
            it_list.append(its)
            x = np.arange(its)
            ax_rob.semilogy(x, error, label=f'dt={dt}')
            ax_rob.legend()
            plt.xlabel('its_num')
            plt.ylabel('log_error')
    fig_rob.savefig(f'error_Robust_C1_{C1}.png')
    ax_tem.semilogx(dts, it_list)
    fig_tem.savefig(f"error_C1_{C1}.png")
    ax.semilogx(dts, it_list, label=f'C1={C1}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
fig.savefig("error_shift.png")







