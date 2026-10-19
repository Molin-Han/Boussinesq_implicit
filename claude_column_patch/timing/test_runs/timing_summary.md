# Timing: star vs single-column (vanka) patch

Minimum over 1 repeat(s) of the PETSc log stage `Second_Run` (one time step, everything compiled and set up once before).  Event times are the max over ranks, inside that stage; PCSetUp and PCApply are nested (outer PC includes the patch work).

| config | ranks | KSP its star -> vanka | time star (s) | time vanka (s) | vanka / star |
|---|---|---|---|---|---|
| nominal | 2 | 13 -> 26 | 23.50 | 28.49 | 1.21 |

| config | ranks | patch | KSP its per Newton step | total KSP | Second_Run (s) | PCSetUp (s) | LU factor (s) | PCApply (s) | MatSolve (s) | time / KSP it (s) |
|---|---|---|---|---|---|---|---|---|---|---|
| nominal | 2 | star | 4/9 | 13 | 23.50 | 2.80 | 0.81 | 22.14 | 9.70 | 1.808 |
| nominal | 2 | vanka | 14/12 | 26 | 28.49 | 2.29 | 0.63 | 26.82 | 8.33 | 1.096 |
