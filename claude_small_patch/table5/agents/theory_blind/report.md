# theory_blind agent report (workflow wf_d02cd389-216, 2026-10-19)

**theory_blind: FS vs MONO for Table 5 (serial runs at nx, nz ≤ 8)**

**Bottom line**
- **Same smoother and coarse grid.** If both multigrids act on the shifted operator A_δ, MONO is algebraically FS. The monolithic column-star patch equals the Schur column-star patch on the velocity error, plus a pressure-residual component that is damped by about 1e-3 per cycle.
- **Measured:** with that setting the iteration counts are equal (8x8: FS 10,11,11 and MONO 10,11,11). With exact inner solves both preconditioners are exactly A_δ^{-1} (11,9,9 each).
- **Cause of the gap:** the Table-5 MONO is not an approximation of A_δ^{-1}. PCMG defaults to `pc_use_amat=TRUE`, and Firedrake's AuxiliaryOperatorPC passes the outer Amat (the unshifted J) to the inner PC. So MONO's finest-level Chebyshev iterates on, and computes residuals with, the unshifted J. Only the patch solves and the coarse levels are shifted.
- **What that buys:** the finest-level smoothing on J removes the column-local (baroclinic) shift error that FS leaves to the outer GMRES. That is the factor of about 2 (5–7 vs 12–13).
- **The paper's explanation is wrong.** `ex_article.tex:3796-3806` says the monolithic patch "resolves more of the pressure-velocity coupling per relaxation sweep". The algebra below and run X1 contradict this.

---

## 1. Operators (from the code)
Unknowns v = (u_xz, u_y, b) and p. Notation:
- C_ij = (w_i·k, b_j), and N² = 1e-4 (`utils.py:19-20, 113-116, 137-139`).
- D_ij = (φ_i, div w_j) (`utils.py:115, 158-161`).

Rotation is off. `--no_rotation` is `store_false` (`LB_time_slice.py:28, 56`), so use_rotation is False unless the flag is passed. The U_mean terms are zero.

The blocks:
- A_vv = [[M_u, 0, -(dt/2)C], [0, M_y, 0], [(dt N²/2)C^T, 0, M_b]]
- J = [[A_vv, -dt D^T], [D, 0]]
- A_δ = J + diag(0, δM_p) (`LB_time_slice.py:204-208`)

Checked numerically: A_vp + dt D^T = 5e-16 (relative), and A_pp = 0.

Define S = A_vv + (dt/δ) D^T M_p^{-1} D. The HDivSchurPC form (`LB_time_slice.py:114-122`) is A_vv + (dt/δ)(div w, div u). This equals S exactly, because the divergence maps extruded RT_2 onto DQ1, which is the pressure space. Numerically, FS with the MatSchurComplement and FS with the auxiliary matrix give the same residual history to 10 digits (`E3_FS_n8` vs `E3_FSnoamat_n8`).

## 2. FS
PCFIELDSPLIT takes every block from Pmat = A_δ:
- `diag_use_amat` and `offdiag_use_amat` default to FALSE (`petsc/src/ksp/pc/impls/fieldsplit/fieldsplit.c:2133-2136`).
- `jac->mat = jac->pmat` (`fieldsplit.c:771`).
- `MatSchurComplementSetSubMatrices` is called with the Pmat blocks (`fieldsplit.c:906`).

The pressure block A00 = δM_p is inverted exactly: bjacobi+ILU(0) of a cell-block-diagonal DQ1 mass matrix is an exact LU. Removing the constant pressure mode does no harm because D^T 1 = 0.

The preconditioner, with S̃^{-1} = one FMG cycle on S:
- P_FS^{-1} = [[I, 0], [-(δM_p)^{-1}D, I]] · diag(S̃^{-1}, (δM_p)^{-1}) · [[I, (dt/δ)D^T M_p^{-1}], [0, I]]
- The output satisfies D x_v + δM_p x_p = r_p exactly.
- Error propagator: I - P_FS^{-1}A_δ = [[E_S, 0], [-(1/δ)M_p^{-1}D E_S, 0]], with E_S = I - S̃^{-1}S. The 1/δ factor amplifies the velocity divergence error into the pressure.

## 3. MONO patch, coarse grid and FMG: the decoupling
**The patch.** The 38Nz+8 star patch contains the 30Nz+8 FS patch plus the pressure of its two columns. That is every pressure dof that D couples to the patch velocity dofs, and M_p is cell-local. So R_i D^T M_p^{-1} D R_i^T = D_i^T M_{p,i}^{-1} D_i, and the local saddle-point solve is the local Schur solve with S_i = R_i S R_i^T followed by local pressure recovery.

**The overlap.** The ASM is BASIC: corrections are summed with no weighting (`firedrake/preconditioners/asm.py:58`). Each cell lies in two stars, so Σ_i R_i^{pT} R_i^p = 2I. The MONO smoother B_M is:
- z_v = B_S (r_v + (dt/δ) D^T M_p^{-1} r_p), where B_S = Σ_i R_i^T S_i^{-1} R_i is exactly FS's level ASM
- z_p = (δM_p)^{-1}(2 r_p - D z_v)
- Equivalently, B_M = P_FS^{-1}[with S̃^{-1} replaced by B_S] + diag(0, (δM_p)^{-1}).

**Change of variables.** Let T = [[I, 0], [(1/δ)M_p^{-1}D, I]], i.e. η = e_p + (1/δ)M_p^{-1}D e_v. This is (δM_p)^{-1} times the shifted continuity residual. Then:
- B_M A_δ = [[B_S S, 0], [(1/δ)M_p^{-1}D(2I - B_S S), 2I]] = T^{-1} diag(B_S S, 2I) T.
- Chebyshev on this operator acts as e_v ← q(B_S S) e_v, which is FS's smoother, and η ← q(2) η.

**Coarse grid.** The spaces are nested, the cells affine and the coefficients constant, so rediscretisation equals Galerkin. The divergence commutes with prolongation: D_h P_v = M_{p,h} P_p M_{p,H}^{-1} D_H. Hence:
- The MONO coarse correction in (e_v, η) is diag(I - P_v S_H^{-1} P_v^T S, I - Π_H), where Π_H is the L2 projection onto coarse DQ1.
- Prolongation maps (e_v, η) to (P_v e_v, P_p η).
- S_H is the same operator as FS's rediscretised coarse S.

**Whole cycle:**
- T (I - P_M,δ^{-1} A_δ) T^{-1} = diag(E_S', E_η), where P_M,δ is MONO with every level on A_δ.
- T (I - P_FS^{-1} A_δ) T^{-1} = diag(E_S, 0).
- E_S' = E_S when the level Chebyshev intervals coincide.
- E_η ≈ q(2)² (I - Π) ≈ 1.3e-3.

**Dense check** (nx=nz=4, n=1344, intervals fixed at [0.2, 2.2]), 2-norms of the blocks of T(I - P^{-1}A_δ)T^{-1}:

| | vv | v-η | η-v | η-η |
|---|---|---|---|---|
| FS | 2.636e-3 (ρ = 1.649e-3) | 3e-12 | 1e-13 | 8e-16 |
| MONO on A_δ | 2.636e-3 | 2e-12 | 3e-13 | 1.291e-3 |

- The vv blocks of FS and MONO-on-A_δ differ by 2.5e-7 relative.
- The η-η norm matches the prediction q(2)² = 1.2905e-3.

**When the two are identical.** With exact inner solves both equal A_δ^{-1}. With one cycle each on A_δ and the same Chebyshev intervals, they differ only by E_η, so counts agree to within ±1–2.

## 4. Chebyshev with PETSc 3.25.3 defaults
Settings (`petsc/src/ksp/ksp/impls/cheby/cheby.c:924-929, 18-38, 800`):
- The estimator is GMRES (CG only if Pmat is flagged SPD), 10 steps on B·Amat with a "noisy" right-hand side.
- emax is the largest real part of the Ritz values; the interval is [0.1, 1.1]·emax.
- Six first-kind steps give max|q| = 0.0478, q(2) = -0.0359 and q'(0) = -9.035.

Estimates at 8x8, levels 1 and 2:
- **FS:** (0.675149, 2.0) and (0.679451, 2.0).
- **MONO on A_δ:** (0.675156, 2.0) and (0.679406, 2.0).
- Both therefore use [0.2, 2.2]. The value 2.0 is exact: it is the eigenvalue of the decoupled u_y mass block (two stars per cell) and of η.
- **MONO as run, finest level:** (0.0922, 1.99335), i.e. [0.1993, 2.1927]. It estimates B_M·J, not B_M·A_δ.
- Fixing [0.2, 2.2] leaves the counts unchanged (E2: 5,5,5).

## 5. What actually differs in the code
**Code path:**
- `firedrake/preconditioners/assembled.py:20, 72`: `pc.setOperators(A, self.P.petscmat)`. The inner PC receives Amat = the outer J.
- `petsc/src/ksp/pc/impls/mg/mg.c:1999`: `pc->useAmat = PETSC_TRUE` in PCCreate_MG.
- `mg.c:935-943`: the finest smoother gets (pc->mat, pc->pmat).
- `mg.c:559-561`: the finest-level residual uses that Amat.
- The ASM is built from Pmat (`asm.py:49`). Levels 0 and 1 use the rediscretised A_δ.
- `params_monolithic['monoshiftpc']` (`LB_time_slice.py:371-381`) does not set `pc_use_amat`, so the default applies.
- The `--inner_rtol` MONO variant uses PCKSP, whose useAmat defaults to FALSE (`pcksp.c:114-116`, `pcksp.c:280-302`), so it is an exact A_δ^{-1}.

**Confirmation from `-ksp_view`** (`runs/MONO_n8/slice3D.txt`):
- The level-2 smoother uses `Mat (firedrake_2_)` (J) as the operator and `Mat (firedrake_2_monoshiftpc_)` (A_δ) to build the preconditioner.
- FS's finest level uses `firedrake_2_fieldsplit_1_` of type schurcomplement, i.e. the shifted S.

**MONO's error propagator on J:**
- I - P_M^{-1}J = q(B_M J)(I - P M1 P^T J) q(B_M J)(I - P M1F P^T J).
- M1 and M1F are the level-1 V-cycle and FMG approximations of A_δ,1^{-1}.
- B_M J = (B_M A_δ)(A_δ^{-1} J).

## 6. How the outer GMRES on J sees each preconditioner
Ordering (v, p):
- **J P_FS^{-1}** = [[S S̃^{-1}, (dt/δ)(S S̃^{-1} - I) D^T M_p^{-1}], [D S̃^{-1}, (dt/δ) D S̃^{-1} D^T M_p^{-1}]].
- **J A_δ^{-1}** (the ideal case S̃ = S) = [[I, 0], [D S^{-1}, (dt/δ) D S^{-1} D^T M_p^{-1}]]. Its pressure eigenvalues are μ = κ/(1+κ), with κ ∈ σ((dt/δ) D A_vv^{-1} D^T M_p^{-1}).

In T-coordinates, with K = [[2S - A_vv, -dt D^T], [(1/δ)M_p^{-1}D, 0]]:
- Ideal: T A_δ^{-1} J T^{-1} = diag(S^{-1}, I) K.
- FS: diag(S̃^{-1}, I) K.
- MONO on A_δ: diag(S̃^{-1}, I - E_η) K.
- MONO's finest smoother: T B_M J T^{-1} = diag(B_S, 2I) K.

**Spectra at 4x4** (`runs/analyse_dense_n4.txt`):
- **Ideal:** 255 pressure eigenvalues differ from 1. 31 are barotropic (vertically uniform), with μ in [4.385e-4, 0.146]; μ1 ≈ (dt/δ)(2π/L)² = 4.4e-4. The other 224 are baroclinic, with μ in [0.33, 0.976]. 159 eigenvalues have |λ-1| > 0.1.
- **FS and MONO on A_δ:** the same outliers; the 12 smallest coincide.
- **MONO as run:** only 18 eigenvalues have |λ-1| > 0.1, and the smallest is 0.0165.

**Mechanism.**
- The patches are whole columns, so on vertical (baroclinic) modes B_M A_δ ≈ 2. That gives B_M J ≈ 2μ ∈ [0.66, 1.95], inside [0.2, 2.2], and those shift errors are removed (factor ≤ 2.1e-3).
- On barotropic modes 2μ ≪ 0.2, so the eigenvalue is 1 - q(2μ)²(1-μ)² ≈ 37μ.

| Barotropic mode | Predicted | Measured |
|---|---|---|
| 1 | 0.0166 | 0.0165 |
| 2 | 0.0649 | 0.0647 |
| 3 | 0.1409 | 0.1403 |
| 4 | 0.2388 | 0.2362 / 0.2382 |

## 7. Predictions and decisive experiments
Prediction: with both preconditioners on A_δ, the counts are equal (±1–2 from E_η). As implemented, MONO needs fewer iterations because it smooths on J.

All runs below were serial, with tmax 200 (counts are per solve). Total compute was about 11.5 minutes. Base commands are as in the task, plus `-firedrake_2_ksp_converged_reason`.

| Run | What it tests | Counts |
|---|---|---|
| X1: MONO + `-firedrake_2_monoshiftpc_pc_use_amat 0` | 8x8 | 10,11,11 (FS 10,11,11; MONO as run 5,5,5) |
| X1 | 4x4 / 4x8 / 8x4 | 9 / 9 / 11, each equal to FS; MONO as run 4 / 4 / 5 |
| X2: `-…_mg_levels_ksp_chebyshev_eigenvalues 0.2,2.2` | prefix `firedrake_2_monoshiftpc_` or `firedrake_2_fieldsplit_1_helmholtzschurpc_` | no change |
| X3: `--inner_rtol 1e-10` | FS vs MONO (outer FGMRES, right) | 11,9,9 and 11,9,9 |
| X3 + `-firedrake_2_monoshiftpc_pc_use_amat 1` | inner solve on J | 1,1,1 |
| X4: MONO `--mono_smoother gmres --outer_ksp fgmres --pc_side right` | default / with `pc_use_amat 0` | 2,2,2 / 11,9,9 (= exact A_δ^{-1}) |
| X5: `--shift` = 1e-3 | FS / MONO on A_δ / MONO as run (8x8) | 22 / 24 / 6 |
| X5: `--shift` = 1e-5 | same | 6 / 6 / 3 |
| X5: `--shift` = 1e-6 | same | 5 / 4 / 3 |

For the production meshes, run X1 and X3 on the Table-5 meshes. Prediction for X1: MONO equals FS (12,13,13,13 left, 11,11,11 right).

## 8. Consequences for the paper
**Choice for the Table-5 MONO line:**
- **(a) Make it the shifted MG.** Add `'pc_use_amat': False` to `monoshiftpc`. It then has the same iteration count as FS with larger patches, so FS wins on time.
- **(b) Keep it.** Describe it as a multigrid on the unshifted J whose patch solves and coarse levels are regularised by δ. The same iterates are then also available with FS-sized 30Nz+8 LU factors, using B_M = P_FS^{-1}[S̃^{-1} → B_S] + diag(0, (δM_p)^{-1}).

**Side findings:**
- **Patch-DOF typo:** `ex_article.tex:3833` lists the (40,80) monolithic patch as 3408 DOFs; it should be 38·80+8 = 3048.
- **Solver name:** the paper says FGMRES, but the Table-5 runs used GMRES (`LB_time_slice.py:330`, `--outer_ksp gmres`).
- **Latent rotation bug (code reading only, not tested):** HDivSchurPC reads `…fieldsplit_1_helmholtzschurpc_use_rotation` (`LB_time_slice.py:106-107`), but the dict sets `fieldsplit_1_use_rotation` (`LB_time_slice.py:317`). With rotation on, the Schur auxiliary operator would drop the Coriolis term. It does not affect Table 5, where rotation is off.

## Files
Everything is in `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_blind/`; no existing file was modified.
- `run.sh` — wrapper for single runs
- `batch_E6E7.sh` — shift and mesh sweeps
- `spectra_driver.py` and `run_dense.sh` — dense P^{-1}, J and A_δ via a KSP monitor
- `pcols_driver.py` — pressure-dof column coordinates
- `analyse_dense.py` — block and spectrum analysis
- `cheb_numbers.py` — Chebyshev polynomial factors
- `runs/<name>/{cmd.txt, log.txt, slice3D.txt}` — per-run command, log and `-ksp_view`
- `runs/analyse_dense_n4.txt`, `runs/cheb_numbers.txt` — analysis output
- `runs/summary_iterations.csv` — all counts and wall times
- `runs/D_*_n4/dense.npz` — 43 MB each; can be deleted
