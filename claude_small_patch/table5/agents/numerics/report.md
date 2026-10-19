# numerics agent report (workflow wf_d02cd389-216, 2026-10-19)

MONO in Table 5 is not the shifted monolithic preconditioner the paper describes. PETSc's PCMG uses the operator matrix Amat by default (`pc->useAmat = PETSC_TRUE`, petsc/src/ksp/pc/impls/mg/mg.c:1999). Firedrake's AuxiliaryOperatorPC passes the outer Amat (the unshifted Jacobian J) to its inner PC (firedrake/preconditioners/assembled.py:72, `pc.setOperators(A, self.P.petscmat)`). So on the finest MG level of MONO, the smoother residual (mg.c:936-939) and the residual sent to the coarse grid (mg.c:1172) both use the unshifted J. Only the star-patch solves (Pmat) and the coarse levels use A_delta.

I checked this at runtime: in MONO the finest smoother's Amat handle is the outer J, and its Pmat is the auxiliary A_delta. FS cannot see J at all: its finest-level Amat is the shifted MatSchurComplement, which equals the auxiliary operator S_aux to 1.7e-23.

With one switch, `pc_use_amat: False` inside `monoshiftpc` (called MONO-P below), MONO behaves like FS:
- the velocity outputs coincide to roundoff;
- with Richardson 0.5 the whole output coincides;
- iteration counts are equal, or MONO-P needs one fewer.

So the user's expectation is correct. The 30%+ gap in Table 5 comes from the PETSc default, not from the pressure treatment.

## Set-up
- **Scripts:** `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/numerics/fs_vs_mono.py`, a standalone copy of LB_time_slice.py:71-212. It imports only utils.py, with bytecode writing off.
- **Problem:** linear Boussinesq slice, first implicit-midpoint step from the initial condition. Rotation off (the default), U_mean = 0, degree 2, bcs and null space as in the script. Length 3e5 m, height 4000 m, dt = 100, shift delta = 1e-4, refinement 2 (horizontal only), nz layers on every level.
- **Meshes:**
  - base nx = nz = 4: fine mesh 16x4 cells, 1344 dofs (u_xz 544, u_y 256, b 288, p 256).
  - (c) only: base nx = nz = 8, fine mesh 32x8 cells, 5248 dofs.
- **FS:** outer KSP with atol 0, rtol 1e-8, max_it 100. Fieldsplit Schur full; split 0 = p (preonly, bjacobi+ilu); split 1 = (u_xz, u_y, b), preonly with HDivSchurPC. Inside that: PCMG `full` cycle, V-cycle, 3 levels; level patch ASMStarPC construct_dim 0 (128 = 30Nz+8 dofs, LU+RCM); coarse LU.
- **MONO:** MonolithicShiftPC → same PCMG on the full A_delta; patch 160 = 38Nz+8 dofs, LU+RCM; coarse LU.
- **Two control solvers:**
  - **FS-P:** FS with `helmholtzschurpc_pc_use_amat false`.
  - **MONO-P:** MONO with `monoshiftpc_pc_use_amat false`.
- **Smoothers on the MG levels (the same in all solvers):**
  - (i) Richardson, scale 0.5, 6 iterations
  - (ii) Chebyshev, 6 iterations, eigenvalues fixed at [0.2, 2.2] (no estimation)
  - (iii) Chebyshev, 6 iterations, PETSc default estimate (the Table-5 setting)
  - (iv) LU on S_aux (FS) or on A_delta (MONO)
  - (v) GMRES(6) smoother with an FGMRES outer solver
- **Exact solves:** LU is PETSc's native LU; Firedrake does not inject MUMPS when `pc_type` is given (solving_utils.py:48-50).
- **Omitted:** `view_patch_sizes` and `ksp_view` (print options only).
- **Compute:** about 4.6 min serial in total; the longest single run was 79 s.

## (a) Auxiliary operator vs. Schur complement
- **Fine level, whole matrix:** ‖S_aux − (A11 − A10 A00⁻¹ A01)‖_F / ‖S‖_F = 1.03e-16.
- **Per block (fine level):**
  - u_xz-u_xz: 1.5e-16, or 2.2e-15 measured against the grad-div part (whose norm is 4.61).
  - u_y-u_y: 1.4e-16.
  - b-b, u_xz-b, b-u_xz: exactly 0.
- **MG levels 0, 1, 2:** the Schur complement of MONO-P's level operator vs FS's level operator differs by 1.07e-16, 1.04e-16 and 1.03e-16 relative.
- **FS's finest-level Amat:** the MatSchurComplement equals the dense S to 1.7e-23.

## Chebyshev intervals in variant (iii)
GMRES estimate, 10 steps, noisy right-hand side; target = [0.1, 1.1] × estimated max.

| solver | level 1: estimate → target | level 2 (finest): estimate → target |
|---|---|---|
| FS (FS-P identical) | (0.667405, 2.0) → [0.2, 2.2] | (0.671272, 2.0) → [0.2, 2.2] |
| MONO-P | (0.667433, 2.0) → [0.2, 2.2] | (0.671231, 2.0) → [0.2, 2.2] |
| MONO (Table 5) | (0.667433, 2.0) → [0.2, 2.2] | (0.124727, 2.00904) → [0.200904, 2.20995] |

- The maximum is exactly 2 in FS: any field supported on the interior dofs of one column is reproduced exactly by both stars that contain that column. In MONO-P the pressure block also gives exactly 2.
- So (iii) is the same as (ii) for FS, FS-P and MONO-P. Only MONO's finest level differs, because there the estimate is taken on M⁻¹J.

## (b) One application of each preconditioner (nx = nz = 4)
- **Norms:** l2 norms of coefficient vectors. "v" means (u_xz, u_y, b). Pressure is compared with its mean removed, because FS's split-0 KSP removes the constant (fieldsplit.c:774-779, kspimpl.h:353-365) and the constant is in ker A.
- **Columns:**
  - dv, dp: difference from FS's output.
  - e_v, e_p: error against x* = A_delta⁻¹ r from a dense LAPACK solve.
  - z = x_p − (δM)⁻¹(r_p − D x_v), relative to ‖x_p‖ (the code's sign convention).
  - eA: error against the unshifted solution x_A = lstsq(J, r); lstsq residual 2.6e-8.
- **Random r:** seed 20261019, bc dofs zeroed.

| variant | solver | dv | dp | e_v | e_p | z |
|---|---|---|---|---|---|---|
| (i) | FS | 0 | 0 | 1.3e-4 | 7.9e-5 | 2.3e-16 |
| (i) | MONO-P | 2.2e-16 | 3.5e-16 | 1.3e-4 | 7.9e-5 | 2.6e-16 |
| (i) | MONO | 0.92 | 1.09 | 0.92 | 1.09 | 0.99 |
| (ii) | FS | 0 | 0 | 8.0e-4 | 8.6e-4 | 2.5e-16 |
| (ii) | MONO-P | 5.0e-16 | 1.28e-4 | 8.0e-4 | 8.6e-4 | 1.28e-4 |
| (ii) | MONO | 0.92 | 1.65 | 0.92 | 1.65 | 0.99 |
| (iii) | MONO-P | 9.4e-14 | 1.28e-4 | 8.0e-4 | 8.6e-4 | 1.28e-4 |
| (iii) | MONO | 0.92 | 1.65 | 0.92 | 1.65 | 0.99 |
| (iv) | MONO | 1.1e-15 | 1.4e-15 | 1.0e-15 | 1.0e-15 | 2.5e-16 |
| (v) | FS | 0 | 0 | 2.3e-10 | 1.8e-10 | 2.3e-16 |
| (v) | MONO-P | 2.4e-16 | 7.6e-14 | 2.3e-10 | 1.8e-10 | 7.6e-14 |
| (v) | MONO | 0.92 | 1.64 | 0.92 | 1.64 | 0.99 |

FS-P matches FS to ≤1.5e-15 in every variant. In (iii), FS gives the same numbers as in (ii).

- **Actual first-step residual:** F(U⁰) with bc rows zeroed. It is identical to the KSP's right-hand side (difference 0.0); ‖r_v‖ = 13.95 and r_p = 0.

| variant | solver | dv | dp | e_v | e_p | z | eA_v | eA_p | unshifted residual after 1 application (v, p) |
|---|---|---|---|---|---|---|---|---|---|
| (i) | FS = MONO-P | 9.6e-17 | 6.1e-16 | 4.0e-6 | 6.9e-6 | 4e-16 | 16.5 | 0.665 | 1.5e-4, 7.48 |
| (i) | MONO | 0.964 | 1.97 | 0.964 | 1.97 | 1.0 | 0.101 | 4.1e-3 | 1.5e-3, 4.6e-2 |
| (ii) | FS = MONO-P | 3.1e-16 | 1.0e-15 | 5.7e-4 | 5.5e-4 | 7e-16 | 16.4 | 0.665 | 5.8e-4, 7.47 |
| (ii) | MONO | 0.968 | 1.98 | 0.968 | 1.98 | 1.0 | 2.33e-2 | 1.27e-3 | 3.0e-3, 1.06e-2 |
| (iii) | FS = MONO-P | 6.8e-14 | 6.7e-14 | 5.7e-4 | 5.5e-4 | 1e-15 | 16.4 | 0.665 | 5.8e-4, 7.47 |
| (iii) | MONO | 0.968 | 1.98 | 0.968 | 1.98 | 1.0 | 2.40e-2 | 1.27e-3 | 3.0e-3, 1.10e-2 |
| (iv) | all | ≤7e-16 | ≤2e-15 | 8e-16 | 1.2e-15 | 3e-16 | 16.5 | 0.665 | 1.5e-14, 7.48 |
| (v) | FS = MONO-P | 1.5e-16 | 5.7e-16 | 2.6e-12 | 4.2e-12 | 4e-16 | 16.5 | 0.665 | 1.6e-10, 7.48 |
| (v) | MONO | 0.970 | 1.99 | 0.970 | 1.99 | 1.0 | 3.3e-6 | 1.9e-7 | 2.7e-6, 1.5e-6 |

In the dv/dp/z columns, "FS = MONO-P" rows give MONO-P's values; FS's own are 0/0/≤4e-16.

**Per field, FS vs MONO-P** (`blockwise_diag.py`): u_xz, u_y and b all agree to ≤2.7e-16 for (i) and (ii), and ≤1.1e-13 for (iii). The pressure differs only in (ii)/(iii) with the random r (1.28e-4). In that case x_p^FS − x_p^MONO-P = −z to 6.5e-12, so the pressure difference is exactly z.

## Mechanism checks (finest level of MONO-P)
- (A_delta M⁻¹ r)_p = 2 r_p holds to 2.5e-15: every column belongs to exactly two stars, and the ASM is BASIC (asm.py:58).
- The ASM velocity output depends only on the reduced residual r_v − A10 A00⁻¹ r_p (difference 1.9e-16).
- One Richardson 0.5 step leaves a pressure residual of 2.5e-15 (zero).
- Chebyshev [0.2, 2.2] × 6 multiplies the pressure residual by −0.035923549746765. The formula q6(2) = T6(−0.8)/T6(1.2) gives −0.035923549746765 (difference 1.2e-13). The PETSc smoother object matches my numpy recurrence to 5e-16.
- Hence z ≠ 0 only when q(2) ≠ 0 and r_p ≠ 0. For the first-step residual r_p = 0, so z is about 1e-15 even with Chebyshev.
- **GMRES caveat:** the agreement for (v) is not structural. GMRES on the full system minimises a different norm, but the pressure part of the preconditioned residual is small here (‖z_v‖ = 572 vs ‖z_p‖ = 0.38). In isolation, one GMRES step differs by 2.2e-9 and six steps by 2.7e-13 (`gmres_diag.py`).

## (c) Outer iteration counts on the unshifted J
GMRES left/right, rtol 1e-8, atol 0; FGMRES (right) for (v).

| variant | FS | FS-P | MONO (Table 5) | MONO-P |
|---|---|---|---|---|
| **4x4 base, 16x4 fine:** (i) | 9/10 | 9/10 | 4/4 | 9/10 |
| (ii) | 9/11 | 9/11 | 4/4 | 9/10 |
| (iii) | 9/11 | 9/11 | 4/4 | 9/10 |
| (iv) | 9/10 | – | 9/10 | – |
| (v) FGMRES | 10 | 10 | 2 | 10 |
| **8x8 base, 32x8 fine:** (i) | 10/11 | 10/11 | 4/3 | 10/11 |
| (ii) | 10/13 | 10/13 | 5/5 | 10/12 |
| (iii) | 10/13 | 10/13 | 5/5 | 10/12 |
| (iv) | 10/11 | – | 10/11 | – |
| (v) FGMRES | 11 | 11 | 2 | 11 |

- FS and MONO-P residual histories with (i) agree to 3e-7 (left) and 4e-9 (right) at every iteration.
- FS(i)'s history matches the exact-LU history (1, 7.9e-3, 7.1e-3, 3.8e-4, …). The count of about 10 is therefore set by the shift δ, not by the MG.
- Left preconditioning stops at different true accuracies: FS true residual 9.2e-7 (4x4) and 2.7e-6 (8x8) with Chebyshev, against 6e-10 for MONO.

## What the numbers support
1. **The auxiliary operator is the exact Schur complement** on every level (about 1e-16).
2. **When the smoother polynomial is the same, FS and the genuinely shifted monolithic MG (MONO-P) give the same velocity** to roundoff in (i), (ii) and (iii). In (iii) this holds because both estimates are exactly 2, so the intervals coincide.
3. **The remaining FS vs MONO-P difference is only in the pressure:** z = −(δM)⁻¹ × (final pressure residual), driven by q(2).
   - It vanishes for Richardson 0.5 (outputs fully identical) and whenever r_p = 0.
   - It changes the counts by at most one iteration, and when it does, MONO-P needs one fewer.
4. **The Table-5 gap is not in the pressure.** Table-5 MONO's velocity differs from FS by 92-97%. Its output approximates the unshifted solution A⁻¹r (eA_v 2.4e-2) rather than A_δ⁻¹r (e_v 0.97). In 4 iterations it beats even the exact shifted inverse (9 iterations). It is a different method: shifted patches and coarse grids, but a finest-level smoother on J.
5. **For the paper:** to compare two shifted preconditioners, rerun the MONO column with `'pc_use_amat': False` in `params_monolithic['monoshiftpc']` (LB_time_slice.py:371-381). The command-line form `-firedrake_2_monoshiftpc_pc_use_amat false` should also work, because that key is not in the dict, but I have not run it through LB_time_slice.py. Otherwise, describe MONO as acting on J at the finest level; its advantage then shows up in true residuals too (right preconditioning: 4-5 vs 10-13 iterations).

Files are in /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/numerics/ (no existing file was modified):
- fs_vs_mono.py
- gmres_diag.py
- blockwise_diag.py
- results_nx4_nz4.json
- results_nx8_nz8_counts_left.json
- results_nx8_nz8_counts_right.json
- run_nx4_nz4.log
- run_nx8_nz8_left.log
- run_nx8_nz8_right.log
- views_nx4_nz4/ (KSP views per solver and level)
