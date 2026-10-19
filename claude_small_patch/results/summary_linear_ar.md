# linear ar test of the paper (delta = C_1/dt): iterations, paper | star re-run | single column

linear: KSP iterations to rtol 1e-8. Folders: data_linear_ar_deg2, data_linear_ar_deg2_star, data_linear_ar_deg2_vanka.

| dt | $L_z = 12$ km | $L_z = 8$ km | $L_z = 4$ km | $L_z = 2$ km | $L_z = 1$ km | $L_z = 0.5$ km | $L_z = 0.25$ km |
|---|---|---|---|---|---|---|---|
| 1 | maxit / maxit / maxit | maxit / maxit / maxit | 80 / 73 / 75 | 41 / 40 / 41 | 18 / 18 / 15 | 6 / 6 / 6 | 5 / 5 / 5 |
| 5 | 78 / 78 / 77 | 52 / 52 / 51 | 28 / 28 / 27 | 12 / 12 / 12 | 7 / 7 / 7 | 5 / 5 / 5 | 5 / 5 / 5 |
| 10 | 52 / 52 / 51 | 37 / 37 / 36 | 18 / 18 / 17 | 9 / 9 / 9 | 6 / 6 / 6 | 5 / 5 / 5 | 4 / 4 / 5 |
| 50 | 19 / 19 / 19 | 14 / 14 / 13 | 8 / 8 / 8 | 6 / 6 / 5 | 5 / 5 / 5 | 4 / 4 / 5 | 4 / 4 / 5 |
| 100 | 13 / 13 / 13 | 10 / 10 / 9 | 6 / 6 / 6 | 5 / 5 / 5 | 4 / 4 / 5 | 4 / 4 / 5 | 4 / 4 / 6 |
| 500 | 8 / 8 / 8 | 6 / 6 / 6 | 5 / 5 / 5 | 4 / 4 / 6 | 4 / 4 / 8 | 4 / 4 / 10 | 4 / 4 / 11 |
| 1000 | 7 / 7 / 7 | 6 / 6 / 6 | 5 / 5 / 5 | 4 / 4 / 5 | 4 / 4 / 7 | 4 / 4 / 9 | 4 / 4 / 11 |
| 5000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 3 / 3 / 4 | 3 / 3 / 5 | 3 / 3 / 6 |
| 10000 | 4 / 4 / 4 | 4 / 4 / 4 | 4 / 4 / 4 | 3 / 3 / 4 | 3 / 3 / 4 | 3 / 3 / 4 | 3 / 3 / 5 |

single column / star over the 61 runs where both converged: median 1.00, mean 1.22, min 0.83, max 2.75; not converged: star 2, column 2; runs still missing: 0

