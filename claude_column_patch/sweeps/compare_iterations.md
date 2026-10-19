# Star patch vs single-column (vanka) patch: iteration counts

Total FGMRES iterations over all Newton steps of the monitored time step, `star / vanka`.  `>100` = a linear solve hit the cap of 100 iterations, `div` = Newton diverged.  Shift $\delta = K\,\delta_{thr}$; $K = 1$ in the dx, dz and ar sweeps.


| sweep | runs compared | median vanka/star | mean | max | not converged (star) | not converged (vanka) |
|---|---|---|---|---|---|---|
| dx | 50 | 1.39 | 1.72 | 5.00 | 0 of 50 | 0 of 50 |
| dz | 50 | 1.63 | 2.18 | 6.50 | 0 of 50 | 0 of 50 |
| ar | 70 | 1.67 | 2.55 | 10.67 | 0 of 70 | 0 of 70 |
| dt | 40 | 1.11 | 1.76 | 5.50 | 10 of 50 | 10 of 50 |


## dx sweep

| dt | N_x (coarse) = 10 | N_x (coarse) = 20 | N_x (coarse) = 40 | N_x (coarse) = 60 | N_x (coarse) = 80 |
|---|---|---|---|---|---|
| 1 | 4 / 19 | 4 / 20 | 4 / 13 | 4 / 13 | 4 / 10 |
| 10 | 6 / 18 | 6 / 18 | 6 / 14 | 6 / 16 | 6 / 14 |
| 50 | 13 / 23 | 13 / 24 | 11 / 19 | 11 / 16 | 11 / 15 |
| 100 | 12 / 25 | 12 / 24 | 10 / 20 | 10 / 18 | 13 / 14 |
| 500 | 16 / 25 | 17 / 29 | 16 / 27 | 16 / 25 | 16 / 22 |
| 1000 | 20 / 29 | 20 / 29 | 20 / 28 | 21 / 27 | 20 / 25 |
| 5000 | 28 / 31 | 29 / 31 | 29 / 30 | 27 / 29 | 27 / 29 |
| 10000 | 28 / 30 | 28 / 29 | 27 / 28 | 25 / 27 | 25 / 27 |
| 50000 | 19 / 21 | 19 / 22 | 19 / 23 | 19 / 23 | 19 / 23 |
| 100000 | 19 / 21 | 18 / 20 | 18 / 21 | 18 / 19 | 18 / 19 |

Ratio vanka/star over the 50 runs where both converged: median 1.39, mean 1.72, max 5.00 (dt = 1, N_x (coarse) = 20: 4 -> 20).  Not converged: star 0, vanka 0 of 50.

Regression against the earlier star sweep (`data_claude_paper_nonlinear_dx`): 50 of 50 common runs identical in total iterations and status.

![dx](compare_dx.png)


## dz sweep

| dt | N_z = 10 | N_z = 20 | N_z = 40 | N_z = 60 | N_z = 80 |
|---|---|---|---|---|---|
| 1 | 6 / 17 | 4 / 18 | 4 / 21 | 4 / 21 | 4 / 26 |
| 10 | 6 / 18 | 6 / 20 | 6 / 18 | 6 / 22 | 6 / 31 |
| 50 | 11 / 20 | 13 / 26 | 11 / 27 | 11 / 32 | 11 / 39 |
| 100 | 12 / 20 | 12 / 24 | 13 / 29 | 13 / 50 | 13 / 40 |
| 500 | 16 / 25 | 17 / 29 | 16 / 38 | 16 / 42 | 16 / 41 |
| 1000 | 20 / 26 | 20 / 29 | 20 / 32 | 21 / 32 | 20 / 36 |
| 5000 | 28 / 30 | 29 / 31 | 29 / 32 | 29 / 33 | 31 / 36 |
| 10000 | 28 / 28 | 28 / 29 | 27 / 30 | 29 / 33 | 29 / 33 |
| 50000 | 19 / 23 | 19 / 22 | 19 / 22 | 19 / 23 | 19 / 23 |
| 100000 | 18 / 20 | 18 / 20 | 18 / 21 | 18 / 21 | 18 / 22 |

Ratio vanka/star over the 50 runs where both converged: median 1.63, mean 2.18, max 6.50 (dt = 1, N_z = 80: 4 -> 26).  Not converged: star 0, vanka 0 of 50.

Regression against the earlier star sweep (`data_claude_paper_nonlinear_dz`): 50 of 50 common runs identical in total iterations and status.

![dz](compare_dz.png)


## ar sweep

| dt | L_z (m) = 12000 | L_z (m) = 8000 | L_z (m) = 4000 | L_z (m) = 2000 | L_z (m) = 1000 | L_z (m) = 500 | L_z (m) = 250 |
|---|---|---|---|---|---|---|---|
| 1 | 8 / 11 | 6 / 13 | 4 / 18 | 4 / 18 | 5 / 28 | 5 / 27 | 3 / 27 |
| 10 | 8 / 13 | 9 / 14 | 6 / 18 | 4 / 29 | 6 / 27 | 5 / 26 | 3 / 20 |
| 50 | 13 / 17 | 10 / 18 | 13 / 24 | 10 / 25 | 6 / 29 | 6 / 29 | 3 / 32 |
| 100 | 13 / 20 | 12 / 21 | 12 / 24 | 11 / 27 | 7 / 33 | 6 / 31 | 5 / 27 |
| 500 | 18 / 25 | 18 / 25 | 17 / 29 | 14 / 33 | 13 / 40 | 11 / 35 | 11 / 39 |
| 1000 | 23 / 27 | 21 / 28 | 20 / 29 | 17 / 33 | 15 / 38 | 13 / 36 | 11 / 35 |
| 5000 | 23 / 24 | 26 / 28 | 29 / 31 | 27 / 35 | 27 / 35 | 20 / 35 | 18 / 37 |
| 10000 | 22 / 23 | 23 / 24 | 28 / 29 | 32 / 34 | 27 / 38 | 27 / 39 | 22 / 38 |
| 50000 | 15 / 19 | 18 / 19 | 19 / 22 | 26 / 28 | 32 / 33 | 37 / 39 | 34 / 44 |
| 100000 | 15 / 19 | 15 / 19 | 18 / 20 | 19 / 23 | 26 / 28 | 32 / 33 | 40 / 44 |

Ratio vanka/star over the 70 runs where both converged: median 1.67, mean 2.55, max 10.67 (dt = 50, L_z (m) = 250: 3 -> 32).  Not converged: star 0, vanka 0 of 70.

Regression against the earlier star sweep (`data_claude_paper_nonlinear_ar`): 70 of 70 common runs identical in total iterations and status.

![ar](compare_ar.png)


## dt sweep

| dt | K = 0.1 | K = 1 | K = 10 | K = 100 | K = 1000 |
|---|---|---|---|---|---|
| 1 | 4 / 12 | 4 / 13 | 6 / 10 | 10 / 10 | 18 / 18 |
| 10 | 4 / 17 | 6 / 18 | 9 / 10 | 16 / 16 | 51 / 51 |
| 50 | 8 / 44 | 13 / 24 | 18 / 21 | 62 / 62 | maxit / maxit |
| 100 | 9 / 30 | 12 / 23 | 21 / 23 | 65 / 61 | maxit / maxit |
| 500 | 11 / 55 | 17 / 29 | 32 / 34 | 106 / 101 | maxit / maxit |
| 1000 | 12 / 53 | 20 / 29 | 40 / 43 | 125 / 132 | maxit / maxit |
| 5000 | 16 / 26 | 29 / 31 | 58 / 58 | maxit / maxit | maxit / maxit |
| 10000 | 16 / 23 | 28 / 29 | 60 / 61 | maxit / maxit | div / div |
| 50000 | 12 / 20 | 19 / 22 | 32 / 32 | 82 / 84 | maxit / maxit |
| 100000 | 12 / 18 | 18 / 20 | 27 / 27 | 55 / 55 | maxit / maxit |

Ratio vanka/star over the 40 runs where both converged: median 1.11, mean 1.76, max 5.50 (dt = 50, K = 0.1: 8 -> 44).  Not converged: star 10, vanka 10 of 50.

Regression against the earlier star sweep (`data_claude_paper_nonlinear_dt`): 50 of 50 common runs identical in total iterations and status.

![dt](compare_dt.png)


## Penalty of the single-column patch against $\rho_1 = m L_z^2/(\pi^2\Delta x^2)$ (all $K = 1$ runs)

| $\rho_1$ | runs | median vanka/star | max vanka/star | median extra iterations | max extra |
|---|---|---|---|---|---|
| [0, 0.01) | 13 | 5.20 | 10.67 | 23 | 29 |
| [0.01, 0.1) | 12 | 3.04 | 7.25 | 16 | 27 |
| [0.1, 1) | 49 | 2.36 | 6.50 | 12 | 37 |
| [1, 10) | 27 | 1.52 | 3.25 | 8 | 16 |
| [10, 100) | 19 | 1.16 | 1.41 | 4 | 11 |
| [100, inf) | 50 | 1.11 | 1.27 | 2 | 4 |

![collapse](collapse_rho1.png)

