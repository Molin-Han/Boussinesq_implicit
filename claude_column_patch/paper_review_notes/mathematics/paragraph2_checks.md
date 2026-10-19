# Referee notes, mathematics lens: second paragraph of Section 5.7 (ex_article.tex lines 3846-3868)

Read-only review, 2026-09-30. Everything below was recomputed or re-derived; file paths are given for each number.
Scripts run: `claude_column_patch/my_checks/check_two_level_constant.py` (venv python) and two short inline
scripts built on `claude_column_patch/theory_kernel/{fe1d,modes}.py` (per-mode two-level minima, discrete
vertical wavenumbers, direct-sum rank).

## (a) "The discrete divergence matrix does not depend on the geometry at all" (line 3846-3848) -- correct

Argument (two sentences). The Raviart-Thomas basis is pushed forward by the contravariant Piola map
u = (det J)^{-1} J u_hat o F^{-1}, for which div u = (det J)^{-1} (div_hat u_hat) o F^{-1} (Piola identity),
while the DG pressure basis is pushed forward by composition, q = q_hat o F^{-1}. Hence
int_K q div u dx = int_Khat q_hat div_hat u_hat dx_hat: every entry of B is a reference-cell integral,
independent of the cell size, of the aspect ratio and even of the cell shape (only the orientation sign enters).

Extruded mesh of `utils.py` / `Nonlinear_Boussinesq_Irk_SC.py` lines 80-83: PeriodicIntervalMesh x uniform
extrusion, re-embedded in R^3 with DG1 vector coordinates (x, 0, z) and a global cell orientation
(`init_cell_orientations(j())`). Every cell is a rectangle dx x dz, J is constant (3 x 2), the immersed Piola map
uses sqrt(det(J^T J)) = dx dz in place of det J and the same cancellation holds; the orientation is the same for
every geometry, so B is identical. `extrude_RT` wraps both tensor-product factors in `HDivElement`, i.e. the
Piola map is applied. Numerical confirmation on the assembled Firedrake operator:
`claude_column_patch/firedrake_algebraic/results_02_kernel_capture.txt`: max|B(Lz=250 m) - B(Lz=4 km)| = 4.2e-26,
max|B(Lz=12 km) - B(Lz=4 km)| = 3.3e-16, max|B(Lx=Lz=1) - B| = 2.2e-16 (max|B| = 0.577), degree 2.

Precision: "only the mass matrix does" -> "only the mass matrices do" (M_u, M_b and M_p all carry the metric;
the kernel of the grad-div term c B^T M_p^{-1} B is ker B whatever M_p is). Optional.

## (b) "spanned by nabla^perp of CG_k(x) x CG_k(z)" (line 3848-3850) -- acceptable, one clarification

With nabla^perp psi = (-psi_z, psi_x): psi_z in CG_k(x) x DG_{k-1}(z) and psi_x in DG_{k-1}(x) x CG_k(z), so
nabla^perp(CG_k x CG_k) is inside the RT space of `extrude_RT(mesh, k)`. The velocity space carries u.n = 0 on top
and bottom, i.e. psi_x = 0 there, i.e. psi constant on the top and constant on the bottom. Hence

  N_h = nabla^perp { psi in CG_k(x) x CG_k(z) : psi constant on top and on bottom }
      = nabla^perp (CG_k(x) x CG_k^0(z))  (+)  span{ nabla^perp z } ,

and nabla^perp z = (-1, 0) is the uniform horizontal flow. Dimension: k N_x (k N_z - 1) + 2 - 1 = k N_x (k N_z - 1) + 1,
which is dim N_h (SVD: `theory_kernel/01_kernel_dimension_counts.txt`, all six cases; Firedrake:
`02_firedrake_crosscheck.txt`, 17 and 81). The uniform flow IS a nabla^perp (of z), only not of a stream function
vanishing on top and bottom (`01_...txt`: "(1,0) contained in curl(Psi_h): False" with Psi_h the zero-trace space).
So "up to the uniform flow" is not needed as long as the stream functions are only required to be constant, not
zero, on top and bottom. As written the sentence omits the boundary condition altogether, and then
nabla^perp(CG_k x CG_k) contains fields with u.n != 0 on top and bottom (dimension k N_x (k N_z + 1) - 1, larger
than dim N_h by 2 k N_x - 2). Also "base mesh" means the 1-D interval mesh elsewhere in the subsection
("closure of a base-mesh cell"), while CG_k(x) x CG_k(z) lives on the x-z slice. Same wording is used at line 1142
of Section 3.4, so this is a consistent simplification; suggested precise form in the findings.

## (c) Counts (k-1)(k N_z - 1) and (2k-1)(k N_z - 1) (line 3851-3853) -- correct

Column patch (closure of base cell j): a divergence-free field of the patch space must have zero u_x on both
lateral face lines, because a lateral-face basis function (vertex hat in x times DG_{k-1}(z)) has a divergence
in the neighbouring column, where the patch has nothing to cancel it. What remains is the kernel of div on one
column with u.n = 0 on all four sides: nabla^perp psi with psi vanishing on the whole column boundary, i.e.
psi in (bubbles of the cell in x, k-1 of them) x CG_k^0(z) (dimension k N_z - 1); the column is simply connected,
so there is no harmonic field. Count (k-1)(k N_z - 1). Star patch of vertex i: two columns, u_x only on the middle
face line, so u.n = 0 on the boundary of the two-column strip: psi in CG_k^0(two cells) x CG_k^0(z), dimension
(2k+1-2)(k N_z - 1) = (2k-1)(k N_z - 1). Checked: `01_kernel_dimension_counts.txt` (k=2, N_z=3: 5 and 15;
k=1: 0 and 2), `02_firedrake_crosscheck.txt` with Firedrake's own patch index sets (same numbers),
`firedrake_algebraic/results_02_kernel_capture.txt` (N_z=4: 7 and 21; N_z=6: 11 and 33; N_z=10: 19 and 57).
N_z = 20, k = 2: 39 and 117.

## (d) CG_2(h) = CG_2(2h) (+) fine bubbles, periodic 1-D mesh coarsened by 2 (line 3855-3856) -- proved

The nodes of CG_2(2h) are the coarse vertices and the coarse midpoints, i.e. exactly the fine vertices. Given
U in CG_2(h) let I_H U in CG_2(2h) interpolate U at the fine vertices; U - I_H U vanishes at every fine vertex, and
a quadratic on a fine cell vanishing at both ends is a multiple of that cell's bubble. So CG_2(h) = CG_2(2h) +
bubbles; the sum is direct because a coarse quadratic that vanishes at all its nodes is zero, and the dimensions
agree: N + N = 2N. Numerically: rank[P | B] = 2N (my inline check: 24 = 12 + 12 at N = 12;
`check_two_level_constant.py`: 32 at N = 16). For k = 1 the complement of CG_1(2h) is spanned by the hats of the
odd vertices (support = two fine cells): no single-cell function exists.

## (e) Coarse kernel + column kernels = whole kernel for k = 2, fails for k = 1; the numbers 81 = 81, 9 < 17 -- correct

Tensorise (d) with CG_2^0(z): CG_2(x;h) x CG_2^0(z) = CG_2(x;2h) x CG_2^0(z) (+) sum_j b_j x CG_2^0(z). Apply
nabla^perp (injective on this space) and add span{nabla^perp z}, which lies in the coarse kernel:
N_h = N_H (+) (+)_j (N_h cap V_j^col). The sum is DIRECT: dim N_H + sum_j dim = [k (N_x/2)(k N_z - 1) + 1] +
(k-1) N_x (k N_z - 1) equals k N_x (k N_z - 1) + 1 exactly when k = 2 (41 + 40 = 81). For k = 1 the column kernels
are empty and dim N_H = 9 < 17 = dim N_h. Requires the hierarchy to coarsen in x only with the same N_z on every
level and nested spaces: `Nonlinear_Boussinesq_Irk_SC.py` line 81-82 (`MeshHierarchy` of the interval, then
`ExtrudedMeshHierarchy(..., layers=[nz]*(refinement+1))`), and `02_firedrake_crosscheck.txt`: the coarse kernel
prolonged by `TransferManager.prolong` has max|div| = 1e-15, Galerkin coarse operator = rediscretised one to
7e-15. The quoted numbers are dim(N_H + sum_j N_h cap V_j) versus dim N_h on N_x = 8, N_z = 3, from
`02_firedrake_crosscheck.txt` (Firedrake's own prolongation and index sets): "[vanka] dim N_H + sum_j = 81,
dim N_h = 81" (k = 2) and "= 9, dim N_h = 17" (k = 1). They are the right numbers; the sentence does not say
which quantity each is (see findings). For k >= 3 the sum still spans but is no longer direct (3.5 N_x(3N_z-1)+1 >
dim N_h); the paper claims k = 2 only in this paragraph, and "k >= 2" only in the closing paragraph, untested.

## (f) Energy of a vertical mode (line 3860-3863) -- correct; sigma and U need one clause of definition

Form of field 1 linearised about rest, f = 0, b eliminated (checked against the assembled fieldsplit_1 matrix to
4e-15 in `firedrake_algebraic/results_01_verify_*`): a(u,w) = (u_x,w_x) + m (u_z,w_z) + c (div u, div w), with
m = 1 + dt^2 N^2 / 4 (paper eq. line 1628). Let sigma_n > 0, g_n in CG_k^0(z), q_n in DG_{k-1}(z) be the singular
system of d/dz with mass-orthonormal g_n, q_n, g_n' = sigma_n q_n (sigma_n^2 = eigenvalues of the discrete
Dirichlet Laplacian in z; plus the barotropic q_0 = const, sigma = 0, which has no u_z). Because every patch is
a whole column and the mesh is a tensor product, modes decouple exactly (`03_mode_decoupling_check.txt`: union
of mode spectra = 2-D spectrum to 1e-14). In mode n: u_x = U(x) q_n(z), u_z = W(x) g_n(z), U in CG_k(x),
W in DG_{k-1}(x); div u = (U' + sigma_n W) q_n, so
a = ||U||^2 + m ||W||^2 + c ||U' + sigma_n W||^2. Divergence free: W = -U'/sigma_n (admissible since
U' in DG_{k-1}(x)); energy = ||U||^2 + (m/sigma_n^2) ||U'||^2 = ||U||^2 + rho dx^2 ||U'||^2 with
rho = m/(sigma_n dx)^2. Correct. The paragraph never says what sigma is (the discrete Dirichlet wavenumber in z)
nor what U is (the horizontal profile of u_x), nor that the operator is the one linearised about rest; the runs
are nonlinear. Suggested wording in the findings.

## (g) The constants 2.0, 12.3, "near 1.3" -- reproduced

`check_two_level_constant.py` (run 2026-09-30): C(0) = 12.2788, C(0.001) = 11.7101, C(0.01) = 8.5045,
C(0.12) = 3.3532, C(1) = 2.1932, C(10) = 2.0199, C(inf) = 2.0000; monotone in rho. Same digits as
`theory_kernel/06_k2_constants.txt` PART 1 (independent code). Definition: C(rho) = sup_U [S(I_H U) + S(U - I_H U)]
/ S(U), S(V) = ||V||^2 + rho h^2 ||V'||^2, I_H = interpolation at the fine vertices. Because the kernel splitting
of (e) is direct, this is the (unique) decomposition constant, and 06 PART 2 shows the two-level additive
lambda_min -> 1/C(rho) as gamma = c/dx^2 -> infinity. At the solver's gamma = 162 (K = 1) the per-mode two-level
lambda_min equals 1/C(rho_n) to within 1 % (my inline computation, all three L_z).

Levels: the same script gives, for the hierarchical splitting with 2, 3, 4, 5 levels, C(0) = 12.28, 15.68, 35.97,
73.57; C(inf) = 2.00, 2.46, 2.67, 2.78. The solver has 3 levels (`--refinement 2`). `04_spectra_parameter_study.csv`,
k = 2, column patch, nominal: two-level lambda_min 0.0820 (1/ = 12.2), three-level 0.0657 (1/ = 15.2);
K = 100: 0.118 (8.5); K = 1000: 0.331 (3.0). So 12.3 is a two-level, delta -> 0 number.

"near 1.3": `04_spectra_parameter_study.csv`, k = 2, star, two-level lambda_min = 0.7878 ... 0.8135 over all
weak-regime cases of the K, L_z, N_x, N_z, dt sweeps, i.e. 1/lambda_min = 1.23 ... 1.27; = 1.00 for
dt >= 5000 s (strong regime); three-level 1.18 ... 1.21. `06_k2_constants.txt` PART 2, star, tau = 60 (thin limit):
0.7937 -> 1.26. "Near 1.3" is a rounding of 1.27; "at most 1.3" is sharper.

## (h) rho_1 = m L_z^2 / (pi^2 dx^2) -- correct

sigma_1 is the smallest discrete Dirichlet wavenumber; for CG_2, N_z = 20, L_z = 4 km: sigma_1 = 7.853985e-4
versus pi/L_z = 7.853982e-4 (ratio 1.000000), so rho_1 = m/(sigma_1 dx)^2 = m L_z^2/(pi^2 dx^2) to seven digits;
nominal dt = 50 s: rho_1 = 0.1225 (the row "0.12" of the study's C(rho) table). Note rho_n DEcreases with n:
the highest mode (n = 39) has rho = 5.1e-5 at the nominal mesh, 4.6e-4 at L_z = 12 km, 2.0e-7 at 250 m.

## (i) Does the paragraph overstate what was proved? -- slightly, in two places

1. "the constant of this splitting ... grows from 2.0 to 12.3" is the constant of one particular splitting
   (unique in the delta -> 0 limit, so it is the two-level decomposition constant), for the additive two-level
   method, for the linearised operator; the three-level number is 15.7, and nothing was derived for the
   multiplicative V-cycle with GMRES smoothing (study record Section 8: "Not proven"). The paragraph should say
   "two-level", "delta -> 0", and, ideally, the three-level value.

2. "thin, weakly stratified domains are the case rho_1 -> 0: the least favourable one". Since rho_n decreases
   with n, the worst mode of the two-level method is the HIGHEST vertical mode, whose rho is ~1e-5 or smaller on
   every mesh of Section 5; the minimum over modes of the two-level lambda_min is 0.0837 (12 km), 0.0820 (4 km),
   0.0818 (250 m), attained at n = 39 in all three cases, and 0.0818 ... 0.0847 over every L_z, N_x, N_z in
   `04_spectra_parameter_study.csv`. The worst-case constant is therefore geometry independent (about 12); what
   rho_1 governs is the constant of the gravest mode, 2.15 (12 km), 3.3 (4 km), 12.0 (250 m), and that mode is
   the one that carries the test data (Remark rem:thin-faster, line 3092). The multigrid-alone experiment with a
   random right-hand side (all modes present) shows the same: 16, 17, 17, 18 iterations for L_z = 12 km ... 250 m,
   study record Section 4.4. The sentence is right in its conclusion but its premise should name the mode.

Everything else in the paragraph checks out.
