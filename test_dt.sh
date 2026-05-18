#!/bin/bash
# run_all.sh
# chmod +x run_all.sh
# to make it runable

# ! Star patch test:

mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 1.0 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 0.1 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 0.01 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 0.001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 0.0001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1.0 --tmax 2.0 --shift 1e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 0.2 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 0.02 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 0.002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 0.0002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 2e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5.0 --tmax 10.0 --shift 2e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 0.1 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 0.01 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 0.001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 0.0001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 1e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10.0 --tmax 20.0 --shift 1e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 0.02 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 0.002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 0.0002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 2e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 2e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 50.0 --tmax 100.0 --shift 2e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 0.01 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 0.001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 0.0001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 100.0 --tmax 200.0 --shift 1e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 0.002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 0.0002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 2e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 2e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 2e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 500.0 --tmax 1000.0 --shift 2e-08 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 0.001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 0.0001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 1e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 1e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 1e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 1000.0 --tmax 2000.0 --shift 1e-08 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 0.0002 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 2e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 2e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 2e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 2e-08 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 5000.0 --tmax 10000.0 --shift 2e-09 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 0.0001 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 1e-05 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 1e-06 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 1e-07 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 1e-08 --maxit 100 --dt_test
mpiexec -n 2 python ../Boussinesq_implicit/LB_time_slice.py --nx 20 --nz 20 --length 300000.0 --height 4000 --dt 10000.0 --tmax 20000.0 --shift 1e-09 --maxit 100 --dt_test
python ../Boussinesq_implicit/plot_shift.py --maxit 100 --dts 1.0 5.0 10.0 50.0 100.0 500.0 1000.0 5000.0 10000.0 --C1_list 1.0 0.1 0.01 0.001 0.0001 1e-05