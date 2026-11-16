levels = ["6"]
imex_dts = [18.75, 37.5, 75, 100]
irk_dts = [14400, 10800, 7200, 3600]
#[37.5, 75, 150, 300, 600, 1200, 2400, 3600]
dts = irk_dts
script = "irk"
pcs = ["al"]
irks = ["RadauIIA"]
stages = [3]
tmax = 86400
ntol = 1.0e-6
ktol = 1.0e-8
williamson=6
ncpus = [16]
gammas = [1.0e3,1.0e4,1.0e5]

warmup = False

import subprocess, os

rows = []

for dt in dts:
    for level in levels:
        for irk in irks:
            for stage in stages:
                for ncpu in ncpus:
                    for pc in pcs:
                        for gamma in gammas:
                            options = {"tmax": tmax,
                                       "ntol": ntol,
                                       "ktol": ktol,
                                       "williamson": 6,
                                       "dt": dt,
                                       "ref_level": level,
                                       "rk_stages": stage,
                                       "pcscheme": pc,
                                       "rk_type": irk,
                                       "gamma": gamma,
                                       }
                            args = []
                            for key, value in options.items():
                                args += ["--"+str(key), str(value)]
                            fname = "irk_data_"+hex(abs(hash(str(options))))
                            try:
                                os.remove(fname)
                            except:
                                pass

                            options["directory"] = fname
                            options["ncpus"] = ncpu
                            os.makedirs(fname)
                            if warmup:
                                args += ["--one_step"]
                                args += ["--show_args"]
                                args += ["--checkpointfile", fname+"/chk.h5"]
                                args += ["--filename", fname+"/data"]
                                args += ["-log_view", ":"+fname+"/log"]
                                args += ["&>", fname+"/out"]
                            print("mpiexec -n "+str(ncpu)+" python irk.py " + " ".join(args))
                            print("grep Main "+fname+"/log &> "+fname+"/stats")
                            print("cat "+fname+"/out >> "+fname+"/stats")

                            rows.append(options)

import pandas as pd
df = pd.DataFrame(rows)
df.to_csv("irks.csv")