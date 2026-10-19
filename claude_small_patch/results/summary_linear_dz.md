# linear dz test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column

linear: KSP iterations to rtol 1e-8. Folders: data_linear_dz_deg2, data_linear_dz_deg2_star, data_linear_dz_deg2_vanka.

| dt | $N_z = 10$ | $N_z = 20$ | $N_z = 40$ | $N_z = 60$ | $N_z = 80$ |
|---|---|---|---|---|---|
| 1 | 50 / 52 / 46 | 80 / 73 / 75 | maxit / maxit / 96 | maxit / maxit / maxit | maxit / maxit / maxit |
| 5 | 27 / 27 / 27 | 28 / 28 / 26 | 28 / 28 / 28 | 28 / 27 / 28 | 28 / 28 / 28 |
| 10 | 18 / 18 / 17 | 18 / 18 / 17 | 18 / 18 / 19 | 18 / 18 / 19 | 18 / 18 / 19 |
| 50 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 |
| 100 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 |
| 500 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 |
| 1000 | 5 / 5 / 5 | 5 / 5 / 5 | 4 / 4 / 5 | 4 / 4 / 5 | 4 / 4 / 5 |
| 5000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 |
| 10000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 3 / 3 / 4 | 3 / 3 / 4 |

single column / star over the 42 runs where both converged: median 1.00, mean 1.03, min 0.88, max 1.33; not converged: star 3, column 2; runs still missing: 0

