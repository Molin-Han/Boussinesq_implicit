#!/bin/bash
# run_all.sh
# chmod +x run_all.sh
# to make it runable

# AR test here:
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 1484 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 545 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 200 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 73 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 27 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 10 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 3 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 1 --dt 100.0 --tmax 200.0 --shift 0.001

# dx test here:
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 85 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 110 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 135 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 160 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 185 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 210 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 235 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 260 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 285 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001

# dz test here:
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 66 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 99 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 147 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 220 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 328 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001

# dt test here:
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1.0 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.1 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.01 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.001 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.0001 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0894427191 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0089442719 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0008944272 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.9443e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.944e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0316227766 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0031622777 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0003162278 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16228e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.1623e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.162e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0028284271 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0002828427 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.82843e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.8284e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.828e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.83e-08 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.0001 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-08 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 0.0003535534 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.53553e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5355e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.536e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.54e-08 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5e-09 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16228e-05 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.1623e-06 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.162e-07 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16e-08 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.2e-09 --dt_test --direct
# mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3e-10 --dt_test --direct





# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.001 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427e-05 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16228e-05 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.8284e-06 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-06 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.536e-07 --dt_test --direct
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16e-08 --dt_test --direct


python ../Boussinesq_implicit/plot_shift.py


