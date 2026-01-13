#!/bin/bash
# run_all.sh

# AR test here:
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 1484 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 545 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 200 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 73 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 27 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 10 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 3 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 1 --dt 100.0 --tmax 200.0 --shift 0.001

# dx test here:
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 85 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 110 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 135 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 160 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 185 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 210 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 235 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 260 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 285 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001

# dz test here:
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 66 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 99 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 147 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 220 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 328 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001

# dt test here:
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1.0
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.089443
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.031623
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.002828
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001
mpiexec -n 16 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.2e-05
