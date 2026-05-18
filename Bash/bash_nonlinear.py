import numpy as np

ncpus = [2] # ! Change to 2 maybe faster.
direct_solver = False # ! if True, use direct solver for Schur complement, else using ASMStar with MG.
maxit = 100

# test = 'dt_test'
# test = 'ar_test'
# test = 'dx_test'
test = 'dz_test'
if test == 'ar_test':
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0]
    heights = [12000, 8000, 4000, 2000, 1000]
    C1_list = [1e-5]
    lengths = [3.0e5]
    nxs = [20]
    nzs = [20]
if test == 'dt_test':
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0]
    C1_list = [0.001, 0.0005, 0.0001, 0.00005, 0.00001, 0.000005, 0.000001]
    heights = [4000]
    lengths = [3.0e5]
    nxs = [20]
    nzs = [20]
if test == 'dx_test':
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0]
    nxs = [10, 20, 40, 60, 80]
    C1_list = [1e-5]
    nzs = [20]
    heights = [4000]
    lengths = [3.0e5]
if test == 'dz_test':
    dts = [1.0, 5.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 5000.0, 10000.0, 50000.0]
    nxs = [20]
    C1_list = [1e-5]
    nzs = [10, 20, 40, 60, 80]
    heights = [4000]
    lengths = [3.0e5]

dts_array = np.array(dts)
tmax_array = dts_array * 2
tmaxs = tmax_array.tolist()

import subprocess, os

rows = []
for nx in nxs:
    for nz in nzs:
        for length in lengths:
            for height in heights:
                for ncpu in ncpus:
                    for dt in dts:
                        for C1 in C1_list:
                            options = {
                                "nx": nx,
                                "nz": nz,
                                "length": length,
                                "height": height,
                                "dt": dt,
                                "tmax": 2 * dt,
                                "shift": np.round(C1 * dt ** (-1), decimals=16),
                            }
                            args = []
                            for key, value in options.items():
                                args += ["--"+str(key), str(value)]
                            fname = "LB_time_data_"+hex(abs(hash(str(options))))
                            try:
                                os.remove(fname)
                            except:
                                pass

                            options["directory"] = fname
                            options["ncpus"] = ncpu
                            # os.makedirs(fname)
                            # if warmup:
                            #     args += ["--one_step"]
                            #     args += ["--show_args"]
                            #     args += ["--checkpointfile", fname+"/chk.h5"]
                            #     args += ["--filename", fname+"/data"]
                            #     args += ["-log_view", ":"+fname+"/log"]
                            #     args += ["&>", fname+"/out"]
                            if direct_solver:
                                print("mpiexec -n "+str(ncpu)+" python ../Boussinesq_implicit/LB_time_slice.py " + " ".join(args)+ ' --maxit' + str(maxit) + ' --' + test + ' --' + 'direct')
                            else:
                                print("mpiexec -n "+str(ncpu)+" python ../Boussinesq_implicit/Nonlinear_Boussinesq_Irk_SC.py " + " ".join(args) + ' --maxit ' + str(maxit) + ' --' + test)
                            # print("grep Main "+fname+"/log &> "+fname+"/stats")
                            # print("cat "+fname+"/out >> "+fname+"/stats")
                            rows.append(options)

# import pandas as pd
# df = pd.DataFrame(rows)
# df.to_csv("irks.csv")

# Plot the errors.

if test == 'dx_test':
    print("python ../Boussinesq_implicit/plot_dx.py " + '--maxit ' + str(maxit) +' --dts ' + " ".join(map(str, dts)) + ' --nxs ' + " ".join(map(str, nxs)) +' --C1 ' + str(C1_list[0]) + ' --length ' + str(lengths[0]))
if test == 'dz_test':
    print("python ../Boussinesq_implicit/plot_dz.py " + '--maxit ' + str(maxit) +' --dts ' + " ".join(map(str, dts)) + ' --nzs ' + " ".join(map(str, nzs)) +' --C1 ' + str(C1_list[0]) + ' --height ' + str(heights[0]))
if test == 'ar_test':
    print("python ../Boussinesq_implicit/plot_ar.py " + '--maxit ' + str(maxit) +' --dts ' + " ".join(map(str, dts)) + ' --heights ' + " ".join(map(str, heights)) +' --C1 ' + str(C1_list[0]) + ' --length ' + str(lengths[0]))
if test == 'dt_test':
    print("python ../Boussinesq_implicit/plot_shift.py " + '--maxit ' + str(maxit) + ' --dts ' + " ".join(map(str, dts)) +' --C1_list ' + " ".join(map(str, C1_list)))