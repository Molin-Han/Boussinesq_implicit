cd "$(dirname "$0")/../../../data/paper_robustness/data_linear_dt/" || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_dt.sh

cd ../data_linear_ar/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_ar.sh

cd ../data_linear_dx/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_dx.sh

cd ../data_linear_dz/ || exit 1
../../../Boussinesq_implicit/sweeps/paper_robustness/test_dz.sh