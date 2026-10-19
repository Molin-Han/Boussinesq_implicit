# linear dx test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column

linear: KSP iterations to rtol 1e-8. Folders: data_linear_dx_deg2, data_linear_dx_deg2_star, data_linear_dx_deg2_vanka.

| dt | $N_x = 10$ | $N_x = 20$ | $N_x = 40$ | $N_x = 60$ | $N_x = 80$ |
|---|---|---|---|---|---|
| 1 | 60 / 62 / 48 | 80 / 73 / 75 | maxit / maxit / maxit | maxit / maxit / maxit | maxit / maxit / maxit |
| 5 | 23 / 23 / 23 | 28 / 28 / 26 | 28 / 28 / 28 | 28 / 28 / 28 | 29 / 29 / 28 |
| 10 | 15 / 15 / 15 | 18 / 18 / 17 | 19 / 19 / 19 | 20 / 20 / 19 | 20 / 20 / 20 |
| 50 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 | 8 / 8 / 8 |
| 100 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 | 6 / 6 / 6 |
| 500 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 |
| 1000 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 | 5 / 5 / 5 |
| 5000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 |
| 10000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 |

single column / star over the 42 runs where both converged: median 1.00, mean 0.99, min 0.77, max 1.03; not converged: star 3, column 3; runs still missing: 0

