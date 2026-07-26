# Flow past a mountain: Boussinesq vertical slice test case

Notes accompanying `Nonlinear_Boussinesq_Mountain.py`.

Reference material:

- C. J. Cotter and J. Shipton, *A compatible finite element discretisation for the
  nonhydrostatic vertical slice equations*, GEM — Int. J. Geomath. **14**:25 (2023),
  [arXiv:2210.07861](https://arxiv.org/abs/2210.07861), section 3.3.
- Reference implementation (compressible Euler):
  [`colinjcotter/sw_implicit/slice_mountain_nh.py`](https://github.com/colinjcotter/sw_implicit/blob/master/slice_mountain_nh.py)
  together with `slice_utils.py`.
- Original test case: Melvin et al. (2010), nonhydrostatic mountain wave.

---

## 1. The test case in the paper (section 3.3, nonhydrostatic case)

A rectangular vertical slice, periodic in `x`, with a rigid lid, whose lower boundary
is raised into a "witch of Agnesi" ridge

$$z_s(x) = \frac{h_m\,a^2}{x^2 + a^2}, \qquad h_m = 1\ \mathrm{m},$$

implemented by the terrain-following coordinate transform

$$(x, z) \;\longmapsto\; \left(x,\; z + z_s(x)\,\frac{H-z}{H}\right).$$

This map leaves the mesh columns vertical (essential for the line / star smoother) and
deforms rectangular cells into trapezia, with the deformation decaying linearly to zero
at the model top.

A uniform horizontal wind is switched on impulsively over a stably stratified
atmosphere at rest. The *only* forcing in the problem is the no-normal-flow condition
`u . n = 0` on the sloping ground, which launches a stationary gravity wave train above
the ridge.

A Newtonian sponge damps the vertical velocity near the lid so that upward-propagating
waves are absorbed rather than reflected:

$$\mu(z) = \begin{cases}
0, & z < z_B, \\[2pt]
\bar\mu \sin^2\!\left(\dfrac{\pi}{2}\dfrac{z - z_B}{H - z_B}\right), & z \ge z_B .
\end{cases}$$

The diagnostic is a contour plot of the vertical velocity at $t = 9000$ s, compared
against Melvin et al. (2010).

---

## 2. Cross-check against `slice_mountain_nh.py`

The repository confirms every parameter, and resolves one contradiction in the paper.

> **Discrepancy.** The text of section 3.3 as rendered states `a = 10000 m` for the
> nonhydrostatic case and `a = 1000 m` for the hydrostatic case. These are **swapped**.
> The code has `a = 1000.` with `U = 10 m/s`, and that is the physically correct
> assignment: "nonhydrostatic" means a *narrow* ridge, and
> $Na/U = 0.01 \times 1000 / 10 = 1$. With `a = 10 km` the 144 km domain would only be
> 14 ridge half-widths wide. **This script uses `a = 1000 m`.**

Everything else agrees between paper and code:

| quantity | value | in the script |
| --- | --- | --- |
| domain length `L` | 144 km | `--length 144.0e3` |
| domain height `H` | 35 km | `--height 35.0e3` |
| columns × layers | 180 × 70 | `--nx 45 --refinement 2` (45·2² = 180), `--nz 70` |
| ridge half width `a` | 1000 m | `--half_width 1000.0` |
| ridge height `h_m` | 1 m | `--h_mount 1.0` |
| mean wind `U` | 10 m/s | `--U_mean 10.0` |
| buoyancy frequency `N` | 0.01 s⁻¹ | `utils.buo_freq()` = 1e-4 = N² |
| timestep | 5 s | `--dt 5.0` |
| final time | 9000 s | `--tmax 9000.0` |
| sponge base `z_B` | H − 10⁴ = 25 km | `--z_sponge 25.0e3` |
| sponge strength | $\bar\mu\,\Delta t = 0.15$ | `--mubar 0.15`, used as `mubar/dt` |
| rotation | none | `--rotation` off by default |

Structural details cross-checked against `slice_utils.py`:

- the sponge enters the momentum equation as `mu * inner(w, k) * inner(u, k) * dx`,
  i.e. Newtonian damping of the **vertical velocity only** — buoyancy/temperature is
  not damped;
- `mu` is defined as `mu_top / dT`, so the paper's $\bar\mu \Delta t$ is the
  timestep-independent quantity;
- boundary conditions are `u . n = 0` on `"top"` and `"bottom"` only (periodic in `x`);
- the time discretisation is the implicit midpoint rule applied to the whole
  nonlinear system (`slice_imr_form`).

Element degrees also match the existing `--degree 2` default. The paper's
next-to-lowest-order compatible spaces

$$V_1 = \mathrm{HDiv}(CG_2 \otimes DG_1) \oplus \mathrm{HDiv}(DG_1 \otimes CG_2), \qquad
V_2 = DG_1 \otimes DG_1, \qquad V_t = DG_1 \otimes CG_2$$

are exactly what `utils.extrude_RT(mesh, k=2)`, `FunctionSpace(mesh, 'DG', 1)` and
`utils.W_theta(mesh, k=2)` produce.

---

## 3. The mathematical problem in the Boussinesq setting

With kinematic pressure $p$ and buoyancy perturbation $b$ about the linear background
$b_{bg} = N^2 z$, the incompressible Boussinesq system solved here is

$$\frac{\partial u}{\partial t} + (u\cdot\nabla)u + f\times u + \mu(z)\,(u\cdot\hat k)\,\hat k
= -\nabla p + b\,\hat k,$$

$$\frac{\partial b}{\partial t} + (u\cdot\nabla)b + N^2 (u\cdot\hat k) = 0,$$

$$\nabla\cdot u = 0,$$

on the periodic slice $0 \le x \le L$, $0 \le z \le H$ with the lower boundary raised by
$z_s$, and $u\cdot n = 0$ on top and bottom. The Coriolis term is switched off for this
test case. The $N^2 (u\cdot\hat k)$ term is already carried by
`utils.Nonlinear_buoyancy_Irk`, and `utils.buo_freq()` already returns
$N^2 = 10^{-4}$, which is the test case value.

**Initial condition:** $u = (U, 0, 0)$, $b = 0$, $p = \mathrm{const}$.

The one structural simplification relative to the compressible version is worth
stating explicitly. In the compressible case the reference script has to bisect on the
Exner boundary value to construct a *discretely* hydrostatic density profile before it
can start. Here no such background solve exists or is needed: with $u$ constant,
$u\cdot\nabla u = 0$, $b = 0$ and $\nabla p = 0$, so the initial state is an **exact
steady solution** of the equations, and the ridge is the sole source of the flow. That
makes the Boussinesq version considerably cleaner as a solver test.

**Regime.** $Na/U = 1$ (fully nonhydrostatic) and $Nh_m/U = 10^{-3} \ll 1$, so the flow
stays in the linear regime where the analytic Queney solution applies. The Boussinesq
and compressible answers should therefore agree, and the expected vertical velocity
scale is

$$w \sim \frac{U h_m}{a} = 10^{-2}\ \mathrm{m/s}.$$

---

## 4. What the script implements

`Nonlinear_Boussinesq_Mountain.py` follows the same skeleton as
`Nonlinear_Boussinesq_Irk_SC.py`: `utils` residual forms, the 3D-embedded slice mesh,
`MeshHierarchy` + `ExtrudedMeshHierarchy`, Irksome `GaussLegendre(1)` (implicit
midpoint), and the `HDivSchurPC(IRKAuxiliaryOperatorPC)` shift preconditioner with the
$p \to p - \delta^{-1}\nabla\cdot u$ elimination, fieldsplit-Schur, a mass solve on the
pressure block and `mg` / `ASMStarPC` on the velocity–buoyancy block.

Four deliberate deviations, all forced by the new physical setting:

### 4.1 Terrain-following mesh hierarchy

`mountain_mesh_hierarchy` replaces `utils.high_dim_mesh_hierarchy`. It embeds each
level in 3D **and** applies the analytic terrain transform on *every* level, so all
levels carry the same orography and retain their column structure. The levels are
geometrically nested only up to the interpolation error of the transform, which is
irrelevant here since $h_m / H \sim 3\times10^{-5}$.

### 4.2 Divergence-free, BC-compatible initial condition

This is the subtlest point. Irksome's default `bc_type="DAE"` imposes both the
algebraic constraint $\nabla\cdot u = 0$ and the Dirichlet data on the **stage values**,
not on $u^{n+1}$ directly. With the implicit midpoint rule the constrained quantity is
$(u^n + u^{n+1})/2$, hence

$$\operatorname{div} u^{n+1} = -\operatorname{div} u^{n},
\qquad u^{n+1}\!\cdot n\big|_{\partial} = -\,u^{n}\!\cdot n\big|_{\partial}.$$

An initial state that is not solenoidal, or that does not respect no-normal-flow on the
slope, is therefore **never cleaned up** — it just flips sign every step. (The
reference `Un -> Unp1` formulation does not have this problem, because there the
boundary condition is imposed on the solution directly.) On a flat mesh the issue is
invisible, because the RT projection of a constant vector is exactly divergence free;
on the deformed mesh it is not.

The script therefore initialises the velocity with a **discrete Leray projection**:
find $(u, p) \in V_1 \times V_2$ with $u\cdot n = 0$ on top and bottom such that

$$\int w\cdot u \,\mathrm{d}x - \int (\nabla\cdot w)\, p \,\mathrm{d}x
= \int w\cdot(U,0,0)\,\mathrm{d}x, \qquad
\int \phi\, \nabla\cdot u \,\mathrm{d}x = 0 .$$

### 4.3 Shift parameter rescaled dimensionally

$\delta$ has units of s/m² — $1/\delta$ is the augmented-Lagrangian parameter
multiplying $B^\top B$ — so the tuned value `1e-4` is a *unit-box* value and does not
transfer to a 144 km × 35 km domain. The default is now

$$\delta = \frac{10^{-4}}{U L} \approx 6.9\times10^{-11},$$

exposed as `--shift_scale` (the dimensionless group $\delta U L$), with a raw override
via `--shift`.

This was **not** cosmetic. With $\delta = 10^{-4}$ the KSP stalled at a residual
reduction rate of 0.99 and Newton failed with `DIVERGED_MAX_IT`; with the rescaled
value the same linear solve converges by `rtol` in 4 iterations.

### 4.4 Solver tolerances and the rotation flag

- `--atol` now defaults to `1e-30`. The mountain is a 1 m perturbation of an exact
  steady state, so absolute residuals are tiny and the KSP has to be driven by the
  relative tolerance; the previous `1e-9` would terminate the linear solves early.
- `--rotation` (default **off**) replaces `--no_rotation`, since the nonhydrostatic
  mountain case is non-rotating.

### 4.5 Diagnostics

The VTK output carries the in-plane velocity, y-velocity, buoyancy, pressure, the
vertical velocity `w` projected into `Vb` (the field plotted in Figure 4 of the paper),
and the cell-wise advective Courant number computed as in the reference script. Max/min
`w` and the max Courant number are printed at every dump.

---

## 5. Cross-checking that it is right

Run at the full paper resolution (180 × 70, `dt = 5 s`, 3 steps):

| check | result |
| --- | --- |
| `\|\|div u₀\|\|_{L²}` after projection | `1.5e-12` |
| normal flow through the terrain, $\int_{\Gamma_b}\lvert u_0\cdot n\rvert$ | `2.9e-15` |
| Newton convergence, step 1 | `1.31 → 4.8e-6 → 3.6e-12` (2 iterations) |
| linear solves | converge by `rtol` in 4–7 iterations |
| **KSP iterations per timestep** | **11.7** |
| `max\|w\|` | `5.6e-3 m/s`, steady across steps |

**Flat-terrain control** (`--h_mount 0`): the initial nonlinear residual is `1.2e-9`,
against `1.31` with the 1 m ridge in place. The discretisation therefore preserves the
uniform-wind / linear-stratification background state to nine digits, and the mountain
signal sits nine orders of magnitude above the numerical floor. (That control run then
aborts, correctly — with nothing left to solve the KSP sits on roundoff. `--h_mount 0`
is only useful for this check.)

**Amplitude check.** The linear scaling gives $U h_m / a = 10^{-2}$ m/s, and Melvin et
al.'s contours for this case span roughly $\pm 5\times10^{-3}$ m/s. The computed
`max|w| = 5.6e-3 m/s` is the right amplitude.

---

## Running it

```bash
# defaults reproduce the paper's nonhydrostatic configuration
python Nonlinear_Boussinesq_Mountain.py

# see all options
python Nonlinear_Boussinesq_Mountain.py --help
```

Compare the `w` field in `Nonlinear_Boussinesq_mountain.pvd` at $t = 9000$ s against
Figure 4 of the paper.

The parameter most worth scanning is `--shift_scale`: it was transferred dimensionally
from the unit-box value rather than tuned at these physical scales.
