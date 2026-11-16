cd "$(dirname "$0")/../../../data/paper_robustness/data_nonlinear_dt_deg2_ew/" || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_nonlinear_dt.sh

cd ../data_nonlinear_ar_deg2_ew/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_nonlinear_ar.sh

cd ../data_nonlinear_dx_deg2_ew/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_nonlinear_dx.sh

cd ../data_nonlinear_dz_deg2_ew/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_nonlinear_dz.sh