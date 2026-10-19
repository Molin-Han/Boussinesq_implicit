# nonlinear dt test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column

nonlinear: total KSP iterations over all Newton steps / Newton steps; f = 0, Eisenstat-Walker on. Folders: data_nonlinear_dt_deg2_ew, data_nonlinear_dt_deg2_ew_star, data_nonlinear_dt_deg2_ew_vanka.

| dt | $C_1 = 10^{-3}$ | $C_1 = 5\times10^{-4}$ | $C_1 = 10^{-4}$ | $C_1 = 5\times10^{-5}$ | $C_1 = 10^{-5}$ | $C_1 = 5\times10^{-6}$ | $C_1 = 10^{-6}$ |
|---|---|---|---|---|---|---|---|
| 1 | missing / div / maxit | missing / maxit / maxit | maxit / maxit / maxit | maxit / maxit / maxit | 91/3 / 90/3 / 90/3 | 63/3 / 64/3 / 64/3 | 33/3 / 32/3 / 32/3 |
| 5 | maxit / maxit / maxit | maxit / maxit / maxit | 82/3 / 88/3 / 88/3 | 59/3 / 62/3 / 62/3 | 29/3 / 30/3 / 30/3 | 21/3 / 22/3 / 22/3 | 18/3 / 13/2 / 13/2 |
| 10 | maxit / maxit / maxit | maxit / maxit / maxit | 87/3 / 105/3 / 105/3 | 63/3 / 71/3 / 71/3 | 29/3 / 33/3 / 33/3 | 22/3 / 25/3 / 25/3 | 13/3 / 8/2 / 8/2 |
| 50 | maxit / maxit / maxit | maxit / maxit / maxit | 72/2 / 73/2 / 73/2 | 57/3 / 54/3 / 54/3 | 24/3 / 22/3 / 26/3 | 14/2 / 14/2 / 17/2 | 13/2 / 10/2 / 20/3 |
| 100 | maxit / maxit / maxit | 129/3 / 133/3 / 131/3 | 34/3 / 34/3 / 35/3 | 29/3 / 30/3 / 30/3 | 17/2 / 14/2 / 19/2 | 15/2 / 12/2 / 24/3 | 11/2 / 8/2 / 26/3 |
| 500 | 33/3 / 28/2 / 29/2 | 26/3 / 27/3 / 28/3 | 21/3 / 19/3 / 23/3 | 16/3 / 12/2 / 24/3 | 9/2 / 10/2 / 28/2 | 9/2 / 9/2 / 25/2 | 8/2 / 7/2 / 29/3 |
| 1000 | 29/3 / 20/2 / 23/2 | 19/3 / 19/3 / 23/3 | 16/3 / 15/3 / 23/3 | 14/3 / 10/2 / 21/2 | 12/3 / 9/2 / 22/2 | 11/3 / 9/3 / 23/2 | 11/3 / 7/2 / 21/2 |
| 5000 | 16/3 / 16/3 / 20/3 | 13/3 / 13/3 / 18/3 | 11/3 / 11/3 / 20/3 | 11/3 / 11/3 / 23/3 | 10/3 / 10/3 / 26/3 | 10/3 / 10/3 / 22/3 | 11/3 / 10/3 / 19/3 |
| 10000 | 14/3 / 14/3 / 18/3 | 13/3 / 13/3 / 15/3 | 12/3 / 11/3 / 19/3 | 11/3 / 10/3 / 19/3 | 11/3 / 10/3 / 21/3 | 11/3 / 10/3 / 20/3 | 11/3 / 10/3 / 26/3 |
| 50000 | 12/3 / 14/3 / 15/3 | 12/3 / 12/3 / 14/2 | 11/3 / 10/3 / 21/3 | 11/3 / 10/3 / 20/3 | 11/3 / 8/2 / 26/3 | 11/3 / 8/2 / 23/3 | 11/3 / 7/2 / 21/3 |

single column / star over the 59 runs where both converged: median 1.29, mean 1.69, min 0.98, max 4.14; not converged: star 11, column 11; runs still missing: 0

