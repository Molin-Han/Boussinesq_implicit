import numpy as np
# nxs = (np.arange(12, 60, 5) * 5).tolist() # ? [60, 85, 110, 135, 160, 185, 210, 235, 260, 285]
# nzs = (np.exp(np.arange(2, 8)* 2/5) * 20) .astype(int).tolist() # ? [44, 66, 99, 147, 220, 328]
# nxs = [60, 80, 100, 150, 200, 250, 300, 350]
nzs = [40, 60, 80, 100, 130, 160, 220, 280, 320]
nxs = [60]
# nzs = [40]
lengths = [3.0e5]
# heights = (np.exp(np.arange(6, -3, -1.0))*10).astype(int).tolist() # ? [4034, 1484, 545, 200, 73, 27, 10, 3, 1]
heights = [4000]
# heights = [12000, 8000, 4000, 2000, 1000, 500, 100, 50]
dts = [1.0, 2.0, 5.0, 8.0, 10.0, 20.0, 50.0, 80.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0] # ! add more time steps
# C1_list = [5.0, 1.0 ,0.1, 0.01, 0.001, 0.0005, 0.0001, 0.00001]
C1_list = [0.001]
# dts = [100.0]
dts_array = np.array(dts)
tmax_array = dts_array * 2
tmaxs = tmax_array.tolist()
# test = 'dt_test'
# test = 'ar_test'
# test = 'dx_test'
test = 'dz_test'

direct_solver = False # ! if True, use direct solver for Schur complement, else using ASMStar with MG.

# shifts = (np.round(dts_array ** (-1.5),decimals=6)).tolist()

ncpus = [2] # ! Change to 2 maybe faster.

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
                                "shift": np.round(C1 * dt ** (-1.5), decimals=16),
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
                                print("mpiexec -n "+str(ncpu)+" python ../Boussinesq_implicit/LB_time_slice.py " + " ".join(args) + ' --' + test + ' --' + 'direct')
                            else:
                                print("mpiexec -n "+str(ncpu)+" python ../Boussinesq_implicit/LB_time_slice.py " + " ".join(args) + ' --' + test)
                            # print("grep Main "+fname+"/log &> "+fname+"/stats")
                            # print("cat "+fname+"/out >> "+fname+"/stats")
                            rows.append(options)

# import pandas as pd
# df = pd.DataFrame(rows)
# df.to_csv("irks.csv")