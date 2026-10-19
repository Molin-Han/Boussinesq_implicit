# linear dt test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column

linear: KSP iterations to rtol 1e-8. Folders: data_linear_dt_deg2, data_linear_dt_deg2_star, data_linear_dt_deg2_vanka.

| dt | $C_1 = 10^{0}$ | $C_1 = 10^{-1}$ | $C_1 = 10^{-2}$ | $C_1 = 10^{-3}$ | $C_1 = 10^{-4}$ | $C_1 = 10^{-5}$ |
|---|---|---|---|---|---|---|
| 1 | maxit / div / div | maxit / maxit / maxit | maxit / maxit / maxit | 80 / 73 / 75 | 26 / 26 / 25 | 12 / 12 / 12 |
| 5 | maxit / maxit / maxit | maxit / maxit / maxit | 64 / 63 / 60 | 28 / 28 / 27 | 12 / 12 / 12 | 6 / 6 / 6 |
| 10 | maxit / maxit / maxit | maxit / maxit / maxit | 48 / 48 / 47 | 18 / 18 / 17 | 9 / 9 / 9 | 5 / 5 / 5 |
| 50 | maxit / maxit / maxit | 48 / 48 / 46 | 17 / 17 / 17 | 8 / 8 / 8 | 5 / 5 / 5 | 4 / 4 / 6 |
| 100 | 90 / 90 / 87 | 33 / 33 / 32 | 13 / 13 / 13 | 6 / 6 / 6 | 4 / 4 / 5 | 4 / 4 / 7 |
| 500 | 48 / 48 / 48 | 19 / 19 / 19 | 8 / 8 / 8 | 5 / 5 / 5 | 4 / 4 / 5 | 3 / 3 / 5 |
| 1000 | 48 / 48 / 48 | 18 / 18 / 18 | 8 / 8 / 8 | 5 / 5 / 5 | 4 / 4 / 5 | 3 / 3 / 4 |
| 5000 | 22 / 22 / 22 | 11 / 11 / 11 | 6 / 6 / 6 | 4 / 4 / 4 | 3 / 3 / 4 | 3 / 3 / 4 |
| 10000 | 15 / 15 / 14 | 8 / 8 / 8 | 5 / 5 / 5 | 4 / 4 / 4 | 3 / 3 / 4 | 3 / 3 / 4 |

single column / star over the 46 runs where both converged: median 1.00, mean 1.09, min 0.93, max 1.75; not converged: star 8, column 8; runs still missing: 0

