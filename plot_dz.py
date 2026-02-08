from firedrake import *
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
from firedrake.output import VTKFile

dts = [1.0, 2.0, 5.0, 8.0, 10.0, 20.0, 50.0, 80.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0]
nzs = [40, 60, 80, 100, 130, 160, 220, 280, 320]
C1 = 0.001
length = 3.0e5
height = 4000

fig, ax = plt.subplots()

for nz in nzs:
    it_list = []
    fig_tem, ax_tem = plt.subplots()
    fig_rob, ax_rob = plt.subplots()
    for dt in dts:
        deltaz = height / nz
        error = np.loadtxt(f'error_dt{dt}_nz{nz}.out')
        its = len(error)
        if its >= 150:
            it_list.append(np.nan)
            # it_list.append(its)
        else:
            it_list.append(its)
        x = np.arange(its)
        ax_rob.semilogy(x, error, label=f'dt={dt}')
        ax_rob.legend()
        plt.xlabel('its_num')
        plt.ylabel('log_error')
    fig_rob.savefig(f'error_Robust_dz_{deltaz}_dt{dt}.png')
    ax_tem.semilogx(dts, it_list)
    fig_tem.savefig(f"error_dz_{deltaz}_dt{dt}.png")
    ax.semilogx(dts, it_list, label=f'dz={deltaz}')
    ax.legend()
    ax.set_xlabel('dt')
    ax.set_ylabel('its')
fig.savefig("error_dz.png")







