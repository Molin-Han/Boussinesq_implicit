import numpy as np
# nxs = (np.arange(12, 60, 5) * 5).tolist() # ? [60, 85, 110, 135, 160, 185, 210, 235, 260, 285]
# nzs = (np.exp(np.arange(2, 8)* 2/5) * 20) .astype(int).tolist() # ? [44, 66, 99, 147, 220, 328]
nxs = [60]
nzs = [44]
lengths = [3.0e5]
# heights = (np.exp(np.arange(6, -3, -1.0))*10).astype(int).tolist() # ? [4034, 1484, 545, 200, 73, 27, 10, 3, 1]
heights = [4034]
dts = [1.0, 5.0, 10.0, 50.0, 100.0, 1000.0]
# dts = [100.0]
dts_array = np.array(dts)
tmax_array = dts_array * 2
tmaxs = tmax_array.tolist()
# shifts = (np.round(dts_array ** (-1.5),decimals=6)).tolist()

ncpus = [16]

import subprocess, os

rows = []
for nx in nxs:
    for nz in nzs:
        for length in lengths:
            for height in heights:
                for ncpu in ncpus:
                    for dt in dts:
                        options = {
                            "nx": nx,
                            "nz": nz,
                            "length": length,
                            "height": height,
                            "dt": dt,
                            "tmax": 2 * dt,
                            "shift": np.round(dt ** (-1.5), decimals=6),
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
                        os.makedirs(fname)
                        # if warmup:
                        #     args += ["--one_step"]
                        #     args += ["--show_args"]
                        #     args += ["--checkpointfile", fname+"/chk.h5"]
                        #     args += ["--filename", fname+"/data"]
                        #     args += ["-log_view", ":"+fname+"/log"]
                        #     args += ["&>", fname+"/out"]
                        print("mpiexec -n "+str(ncpu)+" python ../Boussinesq_implicit/LB_time_slice.py " + " ".join(args))
                        # print("grep Main "+fname+"/log &> "+fname+"/stats")
                        # print("cat "+fname+"/out >> "+fname+"/stats")
                        rows.append(options)

import pandas as pd
df = pd.DataFrame(rows)
df.to_csv("irks.csv")