# theory_refute agent report (workflow wf_d02cd389-216, 2026-10-19)

**theory_refute: findings**

**Main finding.** Claims C2–C4 are correct algebra. They describe a multigrid that runs on A_delta at every level, but the monolithic solver in Table 5 does not do that. PETSc's PCMG uses the outer operator (Amat) by default (`src/petsc/src/ksp/pc/impls/mg/mg.c:1999`, `pc->useAmat = PETSC_TRUE`). Firedrake's AuxiliaryOperatorPC hands its inner PC the outer Amat (`src/firedrake/firedrake/preconditioners/assembled.py:20,72`, `pc.setOperators(A, self.P.petscmat)`). For MONO that Amat is the unshifted Jacobian J.

So on MONO's finest level (`mg.c:936-942`, `mg.c:560`), the Chebyshev smoother, its eigenvalue estimate and the residual passed to the coarse grid all use the unshifted J. A_delta is used only to build the patch solves and the rediscretised coarse levels.

FS goes through the same mechanism, but its outer Amat is PETSc's MatSchurComplement built from the shifted Pmat blocks (`fieldsplit.c:2133-2136`), which equals S. So FS really is a multigrid on S.

The Table 5 run's own `ksp_view` shows this. In `runs/A_mono_np1_nx10_nz10/slice3D.txt:174-201`, the finest-level "linear system matrix" is `Mat (firedrake_2_)` with 549440 nonzeros (J, no pressure-pressure block). The preconditioner matrix is `(firedrake_2_monoshiftpc_)` with 554240 nonzeros (A_delta). In `runs/A_fs_np1_nx10_nz10/slice3D.txt:236-266` the FS finest operator has type `schurcomplement`.

**Experiments:** unmodified `LB_time_slice.py` with the Table 5 flags, serial; iteration counts are for the 3 solves:

| (nx,nz) | FS | MONO as run | MONO + `-firedrake_2_monoshiftpc_pc_use_amat 0` |
|---|---|---|---|
| (4,4) | 9/9/9 | 4/4/4 | 9/9/9 |
| (8,8) | 10/11/11 | 5/5/5 | 10/11/11 |

The same counts come from `probe.py`, which copies the problem and the forms from LB_time_slice.py:
- **(4,4):** FS 9/9, MONO 4/4, MONO with use_amat 0 also 9/9. This holds with estimated and with fixed Chebyshev eigenvalues [0.2, 2.2].
- **(8,8):** FS 10/11, MONO 5/5, MONO with use_amat 0 also 10/11.
- **(8,8), right preconditioning:** FS 13/10, MONO 5/5, MONO with use_amat 0 12/10.

The two inner_rtol (B) runs agree with each other because PCKSP uses Pmat by default (`pcksp.c:114-116`), so both inner solves there are purely shifted.

**C1 (exact inner solves give identical preconditioners and iterations): holds.** In exact arithmetic both preconditioners are A_delta^{-1}. The `fieldsplit_0` solve (bjacobi+ILU) is exact because the DG mass matrix is block diagonal per cell. The constant-pressure null space Firedrake attaches to J, to Jp (`solving_utils.py:354-361`) and to the δM block (`fieldsplit.c:779-781`) only removes a pressure constant that J annihilates and the outer KSP removes anyway (`kspimpl.h:353-365, 414-425`). This matches the B runs: 13/11/11 and 12/10/10 for both solvers.

**C2 (the Schur residual depends only on the velocity error): holds.** It is an identity, and it survives the Dirichlet rows because the boundary rows and columns of B are zero.

**C3 (patch solves give the same velocity correction): holds, for star patches only.** Dense check at (4,4) (N=1344; nv=1088 velocity-like dofs, np=256 pressure dofs):
- 16 patches each, MONO 160 = 38Nz+8 dofs and FS 128 = 30Nz+8 dofs.
- The velocity part of every MONO patch equals the FS patch (16 of 16).
- Every pressure dof is in exactly 2 patches.
- S_aux equals the Schur complement of A_delta to a relative difference of 1.0e-16.

It would fail for ASMVankaPC: its outer-facet velocity dofs reach into cells whose pressure dofs are not in the patch, so the identity R_i^v B^T = B_i^T R_i^p breaks.

**C4 (one MONO cycle gives the same velocity as one FS cycle): the algebra holds, but as a description of the code it fails.**
- The algebra checks out numerically. In ASM_M·A_delta the velocity-pressure block is 1.8e-13 relative, the velocity block equals ASM_S·S to 1.6e-9 (patch conditioning), and the pressure block is 2I to 2.5e-16.
- The full ('full') cycle and coarse-grid argument also hold: refinement is horizontal only, so boundary dofs map to boundary dofs, and restriction is the transpose of prolongation with consistent boundary zeroing (`mg/ufl_utils.py:402-431`, `fmg.c:6-32`).
- With use_amat 0, MONO's and FS's velocity outputs agree to 6e-16 at (4,4) and 6e-13 at (8,8).
- As run, they differ by 84–90% relative, even when the pressure residual r_p is zero.
- Replacing MONO's pressure with FS's back-substituted pressure gives 70/51 iterations at (4,4). MONO's velocity alone is not a shifted-system velocity.

**C5(a) (different Chebyshev intervals because MONO's operator has an extra eigenvalue cluster at 2): fails.**
- ASM_S·S already has eigenvalue exactly 2 for every function supported in one column, because each such function lies in exactly two patches and each reproduces it exactly. At (4,4) that is 960 of 1088 eigenvalues, spectrum real part in [0.668, 2.000].
- MONO's extra pressure eigenvalues at 2 make it 1216 of 1344, with the same interval and the same Euclidean numerical range.
- PETSc's estimate of the largest eigenvalue is exactly 2 for FS and for MONO with use_amat 0, giving identical targets [0.2, 2.2].
- MONO as run estimates 2.006–2.018 on its finest level, but only because it is estimating for ASM_M·J. Fixing both solvers to [0.2, 2.2] leaves the gap at 9 vs 4.

**C5(b) (the pressure output formula): the identity holds, but "the only differences are (a) and (b)" fails.**
- Sharper result: d = (δM)^{-1} × (pressure residual of the shifted system) evolves independently of the velocity error.
- Each Chebyshev smoothing multiplies d by q(2) = T_6(−0.8)/T_6(1.2) = −0.0359; each coarse correction applies (I − Π_H).
- So z = Φ(δM)^{-1} r_p, with Φ ≈ q₂(2)²[(I−Π₁) + q₁(2)⁴(Π₁−Π₀)].
- Measured: |z|/|(δM)^{-1} r_p| = 9.2e-4 at (4,4) and 9.3e-4 at (8,8). z is zero for a velocity-only residual (1e-15) and does not depend on r_v (4e-12).
- Note that Table 5 is left-preconditioned, so the operator GMRES sees is P^{-1}A, not AP^{-1}; the spectra are the same.

**C6 (MONO wins because B^T z cancels the grad-div term): fails.** In a multigrid on A_delta, z depends only on r_p and is about 1e-3 of (δM)^{-1} r_p, so it cannot cancel (1/δ)B^T M^{-1} B e_v. MONO wins because its cycle approximates J^{-1} itself. Dense check at (4,4), first-step right-hand side b (which has b_p = 0):

| | rel. error vs unshifted solution (vel / pres) | ‖Jx − b‖/‖b‖ |
|---|---|---|
| A_delta^{-1} b | 16.4 / 0.665 | 7.48 |
| P_FS^{-1} b (within 5.7e-4 of A_delta^{-1} b) | 16.4 / 0.665 | 7.48 |
| P_MONO^{-1} b, as run | 2.4e-2 / 1.3e-3 | 1.1e-2 |

| eig(P^{-1}J) | min abs | count abs(λ−1) > 0.1 | count abs(λ−1) > 0.01 |
|---|---|---|---|
| FS | 4.4e-4 | 148 | 319 |
| exact A_delta^{-1} | 4.4e-4 | 159 | 255 |
| MONO as run | 1.7e-2 | 18 | 22 |

The smallest exact-shift eigenvalue matches dt·k²/δ for the longest horizontal wavelength: 100·(2π/3e5)²/1e-4 = 4.39e-4. This also explains the 2 iterations with a GMRES(6) level smoother: that smoother minimises the unshifted residual.

**Other things the claims overlook:**
1. **Fix or verify for Table 5:** add `'pc_use_amat': False` to `params_monolithic['monoshiftpc']` (`LB_time_slice.py:371-381`), or pass `-firedrake_2_monoshiftpc_pc_use_amat 0`. It should be rerun at (10,10) through (40,80); I expect FS's counts (12–13).
2. **Wording for the paper:** MONO as run is "the shift used only in the patch solves and the coarse operators; the finest-level smoother acts on the unshifted J". It is not the paper's shifted preconditioner.
3. **Null spaces:** the null spaces attached to the nonsingular A_delta and δM (`assembled.py:80-87`) only project out the pressure constant. They have no effect on the counts.
4. **Noisy right-hand side:** the eigenvalue estimate uses a deterministic hash of the global index with the null space removed (`kspimpl.h:260-283`). The estimate is still pinned at exactly 2 by the eigenvalue-2 cluster, which is why counts are identical across 1, 2 and 16 ranks.
5. **FS 'full' factorisation** (`fieldsplit.c:1262-1322`): the same δM solver is used for both pressure solves. The off-diagonal blocks come from Pmat, but they are identical in J and Jp. The coarse LU uses shift type INBLOCKS (`mg.c:407`), which only acts on zero pivots.
6. **`U_mean` quirk:** `bool(Constant(0.0))` is True, so the `if U_mean:` branches (`utils.py:120-134`, `141-145`) assemble zero-valued advection terms, including interior-facet integrals over dS_v + dS_h. This changes sparsity only, not values, in Table 5. If U_mean were nonzero, velocity advection would be counted twice (upwind block plus centred block).
7. **Rotation option name mismatch:** HDivSchurPC reads `...fieldsplit_1_helmholtzschurpc_use_rotation` (`LB_time_slice.py:107`), but the script sets `fieldsplit_1_use_rotation` (`:317`). This is harmless here because rotation is off; with rotation on, S_aux would differ from the Schur complement of A_delta and C3/C4 would break.

**Files** (all in `/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/theory_refute/`; total compute about 11 minutes, serial, nx ≤ 8):
- Scripts: `probe.py`, `probe3.py`, `probe2.py` (superseded by `probe3.py`), `check_bool.py`
- Logs: `e1_n4.log`, `e2_n4.log`, `e3_n4.log`, `e4_n4.log`, `e5_n8.log`, `e6_n8.log`, `e7a_n8.log`, `e7b_n8.log`, `p3_n4.log`, `p3m_n4.log`
- LB_time_slice.py reruns: `lbrun/lb_*.log`, `lbrun/lb8_*.log` and the matching `lbrun/slice3D*.txt` files
