# code_audit agent report (workflow wf_d02cd389-216, 2026-10-19)

**Main finding.** MONO and FS are the same method in every part except one: which matrix the finest MG level smooths and computes residuals with. In MONO it is the unshifted Jacobian J; in FS it is the shifted Schur complement S_δ. Once MONO is also made to use its shifted auxiliary matrix A_δ there, the two give identical iteration counts in every smoother variant I tried.

How MONO ends up on J:
- `AssembledPC.initialize` hands the inner PC the outer Amat: `pc.setOperators(A, self.P.petscmat)` (firedrake/preconditioners/assembled.py:72; `A, P = pc.getOperators()` at :20).
- `PCCreate_MG` sets `pc->useAmat = PETSC_TRUE` by default (petsc/src/ksp/pc/impls/mg/mg.c:1999).
- So `PCSetUp_MG` gives the finest smoother `(pc->mat, pc->pmat)` (mg.c:935-940), and the finest-level residual operator is that smoother's Amat (mg.c:1168-1173, used in `PCMGMCycle_Private` mg.c:41).
- MONO: the outer PC's Amat is J. The finest smoother iterates on b − J·x with patches from A_δ; coarse levels use A_δ rediscretised.
  - ksp_view, runs/mono_nx8/slice3D.txt and base_mono_nx4/slice3D.txt:174-188: "linear system matrix, followed by the matrix used to construct the preconditioner: (firedrake_2_) … has attached null space / (firedrake_2_monoshiftpc_)".
  - Handle check in the dump: level-2 Amat is the outer J.

How FS ends up on S_δ:
- fieldsplit_1's Amat is the MatSchurComplement (fieldsplit.c:1006), and PCFIELDSPLIT builds every block from the Pmat Jp by default: `jac->mat = jac->pmat` (fieldsplit.c:771), B and C from `pc->pmat` (fieldsplit.c:895-906).
- So FS's finest Amat is S_δ, which is shifted. ksp_view, base_fs_nx4/slice3D.txt:236-240: `type: schurcomplement`.

Proof (serial, fine mesh 32×8, iterations for step 1 / warm-up / official solve):

| Variant | FS | MONO (Table-5 config) | MONO + `-firedrake_2_monoshiftpc_pc_use_amat false` |
|---|---|---|---|
| Chebyshev(6), GMRES left | 10/11/11 | 5/5/5 | 10/11/11 |
| Chebyshev(6), right preconditioning | 13/10/10 | 5/5/5 | 12/10/10 |
| GMRES(6) smoother, FGMRES right | 11/9/9 | 2/2/2 | 11/9/9 |
| Richardson(0.5) | 10/11/11 | 4/3/3 | 10/11/11 |
| Fixed Chebyshev interval [0.2, 2.2] | 10/11/11 | 5/5/5 | 10/11/11 |
| Exact LU inner solve | 10/10/10 | 10/10/10 | – |

On the 16×4 mesh: FS 9/9/9, MONO 4/4/4, MONO without Amat 9/9/9.

- In MONO, `-pc_use_amat false` works from the command line because pcset.c:145-152 reads `-pc_use_amat` after `PCSetType`.
- So one MONO cycle is a multigrid for the unshifted operator that only uses shifted patches and coarse grids. That is why it beats exact A_δ⁻¹ (10–11 iterations here, matching the earlier B runs).
- FS's one cycle approximates A_δ⁻¹.
- The paper's Table 5 MONO numbers (5–7) are therefore not "multigrid on the shifted preconditioner A_δ". With the flag, MONO gives FS's numbers.

**Q1. Is the HDivSchurPC operator exactly the Schur complement of Jp? Yes, with rotation off (as in Table 5).**
- Jacobian of LB_time_slice.py:204-208 using utils.py:113-161:
  - J_uu = ⟨w,u⟩; J_ub = −½dt⟨w·k,b⟩ (bnph); J_up = −dt⟨div w,p⟩ (pressure at n+1, no ½, utils.py:115); J_bu = +½dt N²⟨q,k·u⟩ (utils.py:139); J_bb = ⟨q,b⟩; J_pu = ⟨φ,div u⟩; J_pp = 0.
  - The shift adds +δ⟨p,φ⟩ (LB_time_slice.py:207) using the same `shift` Constant as `appctx["shift"]` (:73, :84). dt is likewise the same Constant (:71, :83).
- The auxiliary form (LB_time_slice.py:114-122) takes lhs() of the same utils forms with uⁿ→One, bⁿ→1 and p := −(1/δ)div u. That gives ⟨w,u⟩ + (dt/δ)⟨div w, div u⟩ − ½dt⟨w·k,b⟩ + ½dt N²⟨q,k·u⟩ + ⟨q,b⟩.
- With A00 = δM_p, A01 = D, A10 = −dt Dᵀ, the term −A10 A00⁻¹ A01 equals (dt/δ)DᵀM_p⁻¹D. This is (div·, div·) exactly when div V ⊂ Q.
- div V ⊂ Q holds:
  - V_2D = HDiv(CG2⊗DG1) ⊕ HDiv(DG1⊗CG2) (utils.py:53-61), so div lands in DG1⊗DG1.
  - The pressure space is DQ1 (the "creating DQ element" warning in the run log).
  - Cells are affine rectangles (constant det J), so the Piola map keeps div in DQ1.
  - div(u_y j) = 0 on the x–z manifold.
- Boundary conditions: `super().form()` returns `ctx.bcs_Jp` of the split problem (base.py:102-113, LB_time_slice.py:124), rebuilt on W1.sub(0) (solving_utils.py:435-442).
- Numerical check, 16×4 mesh:
  - ‖aux − (A11 − A10 A00⁻¹ A01)(Jp)‖_F / ‖S‖_F = 1.03e-16.
  - Identity rows and columns on all 64 BC dofs.
  - The MONO auxiliary matrix equals Jp bit for bit; Jp − J is nonzero only in the pressure block; J_pp = 0.
  - FS's level operator equals the Schur complement of MONO's level operator to about 1e-16 on all 3 levels.
- Rotation is broken when switched on:
  - `--no_rotation` is `store_false` (:28, :56), so passing it turns rotation on.
  - HDivSchurPC reads `…fieldsplit_1_helmholtzschurpc_use_rotation` (:106-107) but the dict sets `…fieldsplit_1_use_rotation` (:317). The auxiliary operator therefore never contains Coriolis.
  - Checked with rotation on: relative error 5.5e-6, max |diff| 19.7, and the difference is skew in the velocity block.
- Side issue: `bool(Constant(0.0))` is True, so `if U_mean:` (utils.py:120, 128, 141) adds zero-valued facet terms to every form. Numerically harmless at U = 0 (it only adds explicit zeros to the sparsity). For U_mean ≠ 0, LB_velocity adds both advection blocks (:120-124 and :127-134): the volume term twice plus two different fluxes.

**Q2. Is fieldsplit_0 an exact solve? Yes, apart from a constant-pressure projection, also in parallel.**
- PETSc first sets the head KSP to GMRES (fieldsplit.c:948); the user's `preonly` then wins at :952. The Schur complement's inner solve is that same KSP (fieldsplit.c:949), and kspupper is too (fieldsplit.c:989).
- The pressure block is 256 rows, 1024 nonzeros (4×4 cell blocks only); ILU(0) reports "fill ratio needed 1." so it is exact.
- Output = (I − 11ᵀ/N)(δM_p)⁻¹x, confirmed to 9.6e-17. The constant is removed because:
  - The solver's nullspace is composed onto W._ises (variational_solver.py:309, solving_utils.py:354-361).
  - The single-field split reuses W._ises[3] (dmhooks.py:376).
  - PCFIELDSPLIT attaches that nullspace to `jac->mat[0]` (fieldsplit.c:774-782).
  - `KSP_PCApply` then calls `KSP_RemoveNullSpace` (kspimpl.h:353-365, 414-425; preonly.c:27).
- This projection is invisible elsewhere: |A10·1| / |A10| = 6e-16, and the outer KSP removes the pressure mean anyway (J carries the nullspace).
- Parallel: pressure dofs belong to the cell's owner and the block is cell-local, so block-Jacobi blocks contain whole cell blocks with no off-process entries. The global fieldsplit_0 KSP still does the mean removal.

**Q3. Are hierarchy, transfers, coarse operators and coarse solves the same apart from pressure? Yes — except for the finest-level Amat above.**
- Levels: 3 (refine level + 1, mg.c:690-696). Horizontal-only coarsening: 4nx, 2nx, nx columns, nz layers on every level.
- `form()` is called once, on the fine level (assembled.py:45). Its UFL form goes into a new context (base.py:149-159); `DMCoarsen` (mg.c:992) → dmhooks.coarsen (dmhooks.py:409-446) → coarsen_snescontext / coarsen_form (mg/ufl_utils.py:85-116, 262-352). That is a symbolic rediscretisation of the same form: coarse arguments, coarse FacetNormal (:56-58), same Constants, coarsened BCs (:126-132).
- Coarse operators are assembled by `compute_operators` (solving_utils.py:527-548) through `setKSPComputeOperators` (dmhooks.py:444-445; DM kept active because Galerkin is "none", mg.c:955, 996-998). Coarse levels have A = P; no nullspace on coarse levels (ufl_utils.py:344-350, the auxiliary context has no `_nullspace`).
- Transfers: TransferManager(use_averaging=True) (solving_utils.py:326).
  - u_xz is a Piola element on a 2D mesh embedded in 3D, so it is non-native (embedded.py:64-66). It goes through DG embedding, native DG transfer and averaged back-projection (embedded.py:263-292); exact for nested spaces.
  - u_y, b and p are native (embedded.py:51-74).
  - Restriction is the exact transpose (embedded.py:311-364). BC dofs are zeroed after prolongation and restriction (ufl_utils.py:402-431).
  - Numerically, MONO's interpolation restricted to (u,b) equals FS's exactly (max diff 0). Field blocks do not couple; pressure rows sum to 1.
- Coarse solve: preonly + LU. Serial: PETSc LU, nd ordering, INBLOCKS shift (mg.c:398-407). Parallel: `lu` on MPIAIJ picks an external package (MUMPS and SuperLU_dist are both available).
- `pc_mg_type full` goes to `PCMGFCycle_Private` (mg.c:635; fmg.c:6-50):
  - restrict the right-hand side to all levels (fmg.c:13-18), x0 = 0 (:24);
  - for each level i < L−1: V-cycle on level i, then overwrite x_{i+1} = I·x_i (:25-31);
  - V-cycle on the finest level (:32).
  - With 3 levels: 3 coarse LU solves, level 1 smoothed 4 × 6 iterations, finest 2 × 6 iterations.

**Q4. Patches.**
- ASMStarPC with construct_dim 0 on the base DMPlex (asm.py:144-188; patches seeded only on owned vertices, :622). The patch is the dofs on the vertical facets above the vertex plus everything inside the two adjacent columns.
- Sizes: FS 30Nz+8 (128 at Nz=4, 248 at Nz=8); MONO 38Nz+8 (160, 312). 4nx patches on the finest level.
- Eliminating the 8Nz pressure dofs from each MONO patch gives exactly the FS patch: error 1.3e-16 over all 24 patches on levels 1 and 2. So the MONO patch is the FS patch plus the pressure dofs of the two columns.
- BC dofs are kept in the patches (only restricted spaces are masked, asm.py:522-530). They appear as identity rows and columns, so they decouple.
- Patch solve: PCASM BASIC (asm.py:58), sorted indices (:66-71), preonly + PETSc LU with RCM and reused ordering, no pivoting.
- LU vs dense LAPACK solve: max relative error 1.3e-15 for both FS and MONO, with patch condition numbers 2.3e8–8.8e8. The local solves are exact.

**Q5. Chebyshev.**
- PETSc version: 3.25.3 development build (release=0), git v3.25.3-375-g9583a52770d, 2026-07-08.
- Defaults (cheby.c:921-929, 375):
  - interval estimated, transform [0, 0.1; 0, 1.1], i.e. targets [0.1·λmax_est, 1.1·λmax_est];
  - 10 estimation steps; noisy right-hand side on; first-kind polynomial.
- Estimator: KSP `<level>esteig_`, GMRES (CG only if Pmat is flagged SPD), rtol 1e-12, max 10 iterations (cheby.c:84-99, 794-801). It reuses the level PC (cheby.c:811) and takes Ritz values from their real parts (cheby.c:18-38).
- The noisy vector is a deterministic hash of the global index, with any Amat nullspace projected out (kspimpl.h:269-284).
- Re-estimated whenever the Amat/Pmat id or state changes (cheby.c:802-856), i.e. at every solve. Each level has its own Chebyshev KSP and estimator.
- Observed on the 16×4 mesh:
  - FS levels 1/2: estimated max 2.0, targets [0.2, 2.2].
  - MONO level 1: [0.2, 2.2]; MONO finest: estimate [0.125, 2.009], targets [0.2009, 2.2099] (because its Amat is J).
  - MONO with Amat off: finest [0.2, 2.2].
  - Fixing [0.2, 2.2] changes nothing in either solver (table above). The Chebyshev interval is not the cause.
- To fix a common interval: `-…mg_levels_ksp_chebyshev_eigenvalues emin,emax` (cheby.c:353-357).
  - This calls `EstEigSet(0,0,0,0)`, which destroys the estimator (cheby.c:74, 109-111); none is recreated (:379, :780). No need to turn esteig off.
  - Do not also pass `-ksp_chebyshev_esteig`: it recreates the estimator (:358-373). Your values still win (:47-60), but the estimate runs for nothing.
  - Constraints: emax > emin, same sign (:69-70).
- ksp_view then prints only "Chebyshev polynomial of first kind / eigenvalue targets used: min, max", with no "estimated via gmres" line and no esteig block (cheby.c:711-745; verified in runs/flag_*_chebfix_nx8).

**Q6. What reaches the MG level smoothers.**
- Prefixes: `firedrake_2_fieldsplit_1_helmholtzschurpc_mg_levels_` (FS) and `firedrake_2_monoshiftpc_mg_levels_` (MONO).
- Identical options from the dict: ksp_type chebyshev, ksp_max_it 6, pc_type python firedrake.ASMStarPC, pc_star_construct_dim 0, pc_star_sub_sub_pc_type lu, view_patch_sizes, plus RCM and reuse_ordering. Also mg_coarse preonly/lu, pc_mg_type full, cycle v, auxiliary mat_type aij (FS by default, assembled.py:42; MONO explicit).
- Identical implicit PCMG settings:
  - norm type NONE plus KSPConvergedSkip (mg.c:412-413; ksp_view says "not checking for convergence");
  - nonzero initial guess (mg.c:1158-1162);
  - the same KSP for down- and up-smoothing (mg.c:436), so 6 iterations pre and 6 post;
  - left preconditioning.
- Differences beyond the pressure field:
  1. **Finest Amat:** J (unshifted) in MONO vs S_δ as a MatSchurComplement in FS. This is the cause.
  2. **Nullspace:** MONO's finest Amat (J) and its auxiliary matrix (copied from Jp, assembled.py:58-59, 80-88) carry the constant-pressure nullspace. The Chebyshev smoother there therefore removes the pressure mean after every PC apply and projects the noisy vector. FS's level matrices have no nullspace (fieldsplit.c:913-914); FS removes the mean only inside A00⁻¹.
  3. **Structure:** FS runs the cycle on r_ub − A10 A00⁻¹ r_p and then recovers the pressure exactly; MONO runs one cycle on the full system.
  4. **Cost:** FS's finest-level MatMult also runs a fieldsplit_0 solve.
  5. **Unused dict key:** FS has 'fieldsplit_1_use_rotation'.
- Outer solver: GMRES left in both, but the convergence test is on the residual preconditioned by each solver's own preconditioner, so the two measure different norms.

**Q7. Flags to add to LB_time_slice.py.** All tested in audit/lb_audit_copy.py; default behaviour unchanged.

Insert after line 48:
```python
parser.add_argument('--inner_lu', action='store_true', help='Diagnostic: exact LU inner solves in both solvers (FS: LU of the assembled Schur auxiliary operator S_delta instead of the MG cycle; MONO: LU of A_delta instead of the MG cycle).')
parser.add_argument('--fs_smoother', type=str, default=None, choices=['chebyshev', 'richardson', 'gmres'], help='Diagnostic: KSP around the patch smoother on the MG levels of the fieldsplit solver, with fieldsplit_1 kept preonly (ONE MG cycle per outer iteration). Default: chebyshev with --richardson, gmres (and fieldsplit_1 fgmres) without.')
parser.add_argument('--richardson_scale', type=float, default=None, help='Diagnostic: damping of a Richardson MG level smoother, both solvers (PETSc default 1.0).')
parser.add_argument('--cheb_interval', type=float, nargs=2, default=None, metavar=('EMIN', 'EMAX'), help='Diagnostic: fixed Chebyshev interval on every MG level of both solvers (-ksp_chebyshev_eigenvalues; switches the GMRES eigenvalue estimate off).')
parser.add_argument('--mg_use_pmat', action='store_true', help='Diagnostic: finest-level smoother and residual of the MG cycle on the auxiliary (shifted) operator instead of the outer Amat (PCMG default useAmat=TRUE: unshifted J for MONO, the shifted Schur complement for FS).')
```

Insert immediately before line 294 (`params_schur = {`):
```python
if args.inner_lu:   # (i) exact LU of the assembled S_delta instead of the MG cycle
    helmholtz_schur_pc_params = {'pc_type': 'lu', 'pc_factor_mat_solver_type': 'mumps'}
elif not args.direct:
    if args.fs_smoother is not None:          # (ii)/(iv) level KSP of the FS cycle
        helmholtz_schur_pc_params['mg_levels_ksp_type'] = args.fs_smoother
    if args.richardson_scale is not None:     # (ii)
        helmholtz_schur_pc_params['mg_levels_ksp_richardson_scale'] = args.richardson_scale
    if args.cheb_interval is not None:        # (iii)
        helmholtz_schur_pc_params['mg_levels_ksp_chebyshev_eigenvalues'] = '%r,%r' % tuple(args.cheb_interval)
    if args.mg_use_pmat:                      # (v)
        helmholtz_schur_pc_params['pc_use_amat'] = False
```

Insert after line 331 (before the `--inner_rtol` block, which still wins):
```python
if args.fs_smoother is not None:
    params_schur['fieldsplit_1_ksp_type'] = 'preonly'
```

Insert after line 358:
```python
if args.richardson_scale is not None:         # (ii)
    mg_levels_mono['ksp_richardson_scale'] = args.richardson_scale
if args.cheb_interval is not None:            # (iii)
    mg_levels_mono['ksp_chebyshev_eigenvalues'] = '%r,%r' % tuple(args.cheb_interval)
```

Insert after line 391 (before `if args.pc_side`):
```python
if args.inner_lu:   # (i) exact LU of A_delta instead of the MG cycle
    params_monolithic['monoshiftpc'] = {'mat_type': 'aij', 'pc_type': 'lu', 'pc_factor_mat_solver_type': 'mumps'}
elif args.mg_use_pmat and args.inner_rtol is None:   # (v)
    params_monolithic['monoshiftpc']['pc_use_amat'] = False
```

How to run each case (MONO = `--monolithic --outer_ksp gmres --pc_side left`):
- **(i) Exact LU:** FS `--richardson --inner_lu` (same as the existing `--direct`); MONO `… --mono_smoother chebyshev --inner_lu`.
- **(ii) Richardson of scale s:** FS `--richardson --fs_smoother richardson --richardson_scale s`; MONO `… --mono_smoother richardson --richardson_scale s`.
- **(iii) Fixed Chebyshev interval:** add `--cheb_interval EMIN EMAX` to either solver.
- **(iv) FS with GMRES smoother, one cycle:** `--fs_smoother gmres` without `--richardson` (outer becomes FGMRES), or `--richardson --fs_smoother gmres --outer_ksp fgmres`. The MONO counterpart is `--monolithic --mono_smoother gmres`.
- **(v) Make MONO the same method as FS:** add `--mg_use_pmat`.

Command-line equivalents that need no code change (these keys are not in the dict, so they are not overridden):
- `-firedrake_2_monoshiftpc_pc_use_amat false` and `-firedrake_2_fieldsplit_1_helmholtzschurpc_pc_use_amat false` (both tested).
- `-firedrake_2_{monoshiftpc|fieldsplit_1_helmholtzschurpc}_mg_levels_ksp_chebyshev_eigenvalues EMIN,EMAX` (tested).
- `-firedrake_2_monoshiftpc_mg_levels_ksp_richardson_scale s` (with `--mono_smoother richardson`).

Anything that changes a key already in the dict (FS level ksp_type, fieldsplit_1_ksp_type, monoshiftpc pc_type) needs the flags: auto-prefixed dict keys are re-inserted and then deleted around every solve (petsctools/options.py:418-423, 526-537). The `_2_` in the prefix is a global counter: the two `project()` calls at :178-179 take `_0_` and `_1_`, so adding any solve before line 421 shifts it.

**Experiment parameters (all runs).** Linear problem (LB_time_slice.py: linearised Boussinesq vertical slice, implicit midpoint). Serial, no mpiexec. `--length 300000 --height 4000 --dt 100 --tmax 200` (3 solves: step 1, warm-up, official), `--shift 1e-4 --rtol 1e-8 --maxit 100 --reordering --timing`, degree 2, refinement 2 (3 MG levels, fine mesh 4nx × nz, i.e. 32×8 for `--nx 8 --nz 8` or 16×4 for `--nx 4 --nz 4`). Star patch (ASMStarPC construct_dim 0, LU+RCM), MG full V-cycle, 6 smoother iterations, coarse LU. Outer GMRES left unless stated. Total compute was about 8 minutes.

Everything is in /home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/code_audit:
- run_audit.sh — runs one configuration; outputs in runs/<tag>/{cmd.sh,run.log,slice3D.txt}
- base_fs_nx4/, base_mono_nx4/, mono_noamat_nx4/, fs_direct_nx4/ — the first 16×4 runs (ksp_view in slice3D.txt)
- audit/make_copy.py — builds audit/lb_audit_copy.py from LB_time_slice.py plus the Q7 flags and an operator dump
- audit/analyse_dumps.py and its output audit/analyse_nx4.txt — the numerical checks for Q1–Q4
