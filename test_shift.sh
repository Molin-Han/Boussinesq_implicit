#!/bin/bash
# run_all.sh
# chmod +x run_all.sh
# to make it runable

# ! Star patch test:
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 5.0 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1.0 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.1 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.01 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.001 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.0005 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.0001 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 1.7677669529663689 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.3535533905932738 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0353553390593274 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0035355339059327 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0003535533905933 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0001767766952966 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 3.53553390593e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 3.5355339059e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.447213595499958 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0894427190999916 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0089442719099992 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0008944271909999 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 4.472135955e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.2209708691207961 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0441941738241592 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0044194173824159 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0004419417382416 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.41941738242e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 2.20970869121e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.4194173824e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.419417382e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.158113883008419 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0316227766016838 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0031622776601684 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0003162277660168 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16227766017e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 1.58113883008e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.1622776602e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16227766e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0559016994374947 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0111803398874989 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0011180339887499 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.000111803398875 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.11803398875e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 5.5901699437e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.1180339887e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.118033989e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.014142135623731 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0028284271247462 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0002828427124746 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.82842712475e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.8284271247e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 1.4142135624e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.828427125e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.82842712e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0069877124296868 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0013975424859374 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0001397542485937 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.39754248594e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.3975424859e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 6.98771243e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.397542486e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.39754249e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.005 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.0001 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 5e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 0.0017677669529664 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 0.0003535533905933 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.53553390593e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5355339059e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.535533906e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 1.767766953e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.53553391e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5355339e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 0.0004472135955 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 4.47213595e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.9442719e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.944272e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 0.0001581138830084 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766017e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.1622776602e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 1.58113883e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.1622777e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.162278e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 5.59016994375e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.11803398875e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.1180339887e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118033989e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.11803399e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 5.5901699e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118034e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118034e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 1.41421356237e-05 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.8284271247e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.828427125e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.82842712e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.8284271e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 1.4142136e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.828427e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.82843e-11 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 5e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 5e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-11 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 1.767766953e-06 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.535533906e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.53553391e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.5355339e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.535534e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 1.767767e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.53553e-11 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.5355e-12 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 4.472135955e-07 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.94427191e-08 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.9442719e-09 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.944272e-10 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.94427e-11 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 4.47214e-11 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.9443e-12 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.944e-13 --dt_test



# ! Direct test:
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 5.0 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1.0 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.1 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.01 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.001 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.0005 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 0.0001 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1.0 --tmax 2.0 --shift 1e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 1.7677669529663689 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.3535533905932738 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0353553390593274 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0035355339059327 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0003535533905933 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 0.0001767766952966 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 3.53553390593e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2.0 --tmax 4.0 --shift 3.5355339059e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.447213595499958 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0894427190999916 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0089442719099992 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 0.0008944271909999 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 4.472135955e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5.0 --tmax 10.0 --shift 8.94427191e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.2209708691207961 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0441941738241592 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0044194173824159 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 0.0004419417382416 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.41941738242e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 2.20970869121e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.4194173824e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 8.0 --tmax 16.0 --shift 4.419417382e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.158113883008419 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0316227766016838 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0031622776601684 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 0.0003162277660168 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16227766017e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 1.58113883008e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.1622776602e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10.0 --tmax 20.0 --shift 3.16227766e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0559016994374947 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0111803398874989 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.0011180339887499 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 0.000111803398875 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.11803398875e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 5.5901699437e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.1180339887e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20.0 --tmax 40.0 --shift 1.118033989e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.014142135623731 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0028284271247462 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 0.0002828427124746 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.82842712475e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.8284271247e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 1.4142135624e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.828427125e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50.0 --tmax 100.0 --shift 2.82842712e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0069877124296868 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0013975424859374 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 0.0001397542485937 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.39754248594e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.3975424859e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 6.98771243e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.397542486e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 80.0 --tmax 160.0 --shift 1.39754249e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.005 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.001 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 0.0001 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 5e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 100.0 --tmax 200.0 --shift 1e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 0.0017677669529664 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 0.0003535533905933 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.53553390593e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5355339059e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.535533906e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 1.767766953e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.53553391e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 200.0 --tmax 400.0 --shift 3.5355339e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 0.0004472135955 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.94427191e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 4.47213595e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.9442719e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 500.0 --tmax 1000.0 --shift 8.944272e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 0.0001581138830084 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766017e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.1622776602e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.16227766e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 1.58113883e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.1622777e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 1000.0 --tmax 2000.0 --shift 3.162278e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 5.59016994375e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.11803398875e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.1180339887e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118033989e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.11803399e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 5.5901699e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118034e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 2000.0 --tmax 4000.0 --shift 1.118034e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 1.41421356237e-05 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.8284271247e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.828427125e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.82842712e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.8284271e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 1.4142136e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.828427e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 5000.0 --tmax 10000.0 --shift 2.82843e-11 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 5e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 5e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 10000.0 --tmax 20000.0 --shift 1e-11 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 1.767766953e-06 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.535533906e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.53553391e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.5355339e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.535534e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 1.767767e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.53553e-11 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 20000.0 --tmax 40000.0 --shift 3.5355e-12 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 4.472135955e-07 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.94427191e-08 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.9442719e-09 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.944272e-10 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.94427e-11 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 4.47214e-11 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.9443e-12 --dt_test --direct
# mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 60 --nz 44 --length 300000.0 --height 4034 --dt 50000.0 --tmax 100000.0 --shift 8.944e-13 --dt_test --direct


python ../Boussinesq_implicit/plot_shift.py


