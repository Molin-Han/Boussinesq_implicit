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

> **These checks were not sufficient.** They test that the *solver* converges to the
> solution of the discrete system. They say nothing about whether the discrete system
> is the right one. Both bugs in section 6 leave every one of these numbers looking
> healthy. The amplitude check that used to be in this section (`max|w| = 5.6e-3 m/s`,
> "the right amplitude") was reading a contaminated field: it is the correct order of
> magnitude by coincidence, because the spurious forcing happens to scale with the
> same $U h_m/a$. **A converged residual is not evidence of a correct discretisation** —
> the checks that actually caught these were the operator-consistency test and the
> mesh-continuity test in section 6.

**Flat-terrain control** (`--h_mount 0`): the initial nonlinear residual is `1.2e-9`,
against `1.31` with the 1 m ridge in place. This control is still meaningful: it shows
the discretisation preserves the uniform-wind / linear-stratification background state
to nine digits. Note it does *not* catch either bug — with `h_mount = 0` the mesh is
undeformed (no tear) and $|u|^2$ is constant (the spurious advection force cancels).

---

## 6. Two bugs found after the first version, and their fixes

The symptom was a velocity field with visible cell-to-cell jumps, while the solver
converged perfectly. That combination points away from the solver and towards the
discretisation. Two independent causes were found; both are verified in isolation
below.

### 6.1 The nonlinear advection operator was transposed — `utils.py`

```python
# utils.Nonlinear_velocity_Irk, as written
eqn -= inner(div(outer(u, w)), u) * dx      # WRONG
# slice_utils.u_tendency (Cotter), and the fix
eqn -= inner(div(outer(w, u)), u) * dx      # right
```

UFL contracts `div` over the **last** index, so `outer` is not symmetric here:

$$\operatorname{div}(w \otimes u)_i = \partial_j (w_i u_j) = (u\cdot\nabla)w_i + w_i \nabla\!\cdot\! u,$$

$$\operatorname{div}(u \otimes w)_i = \partial_j (u_i w_j) = (w\cdot\nabla)u_i + u_i \nabla\!\cdot\! w .$$

Only the first is the integration by parts of $\int w\cdot(u\cdot\nabla)u$.

**Verification.** On a smooth, continuous, divergence-free field with $u\cdot n = 0$ on
top and bottom and periodic in $x$ — where every facet and boundary term vanishes, so
the weak and strong forms must agree exactly:

| form | relative difference from $\int w\cdot(u\cdot\nabla)u$ |
| --- | --- |
| `-∫ div(outer(w,u))·u` (Cotter, the fix) | `2.6e-14` ✅ |
| `-∫ div(outer(u,w))·u` (utils.py) | `7.4e+00` ❌ |

**Why it produces jumps at cell interfaces.** Writing $K = \lvert u\rvert^2/2$, the
erroneous term expands to

$$-\!\int\! (w\cdot\nabla)K - \!\int\! K\,\nabla\!\cdot\! w
\;=\; -\!\int\! \nabla\!\cdot\!(wK) \;-\; \int\! K\,\nabla\!\cdot\! w .$$

The second piece is pressure-shaped and is harmlessly absorbed by the pressure Lagrange
multiplier. The first, summed over cells, is $-\sum_K \oint_{\partial K} K\, w\cdot n$,
and since $w \in H(\mathrm{div})$ has continuous normal component this collapses to a
spurious **interfacial force $\propto [\![K]\!]\, w\cdot n$ on every facet** — a forcing
that lives exactly on cell boundaries. With $U = 10$ m/s, $K \approx 50$ m²/s², and the
jump in $K$ across facets comes from the tangential velocity jumps that an
$H(\mathrm{div})$ space legitimately has. The result swamps the $5\times10^{-3}$ m/s
mountain signal.

**Why it never showed up before.** In the earlier tests `U_mean = 0` and
$\lvert u\rvert \sim 10^{-3}$ m/s, so $K \sim 10^{-6}$ and the whole term — right or
wrong — is negligible; those runs were effectively linear. The mountain case is the
first one with a real mean flow, and there the advection term is the term that carries
the entire solution.

**Scope.** The same transposition was present in three places in `utils.py`, all now
fixed:

- `Nonlinear_velocity_Irk` — used by this script, `Nonlinear_Boussinesq_Irk_SC.py`,
  `Nonlinear_Irk_Boussinesq_direct.py`;
- `Nonlinear_velocity` — used by `Nonlinear_Boussinesq_slice.py`;
- `LB_velocity`, in the `U_mean` mean-flow branch (twice — the upwind block and the
  centred-flux block). Inactive at the default `U_mean = 0`.

Anything previously run with a **nonzero mean flow or a finite-amplitude velocity** is
affected. The small-amplitude, zero-mean-flow convergence and scaling studies are not:
there the term is quadratically negligible.

### 6.2 The terrain-following mesh was torn open — mesh construction

`utils.high_dim_mesh_hierarchy` builds the embedded mesh by interpolating the new
coordinates into `VectorFunctionSpace(m, "DG", 1)`. Firedrake's default `DG` variant is
`spectral`, whose degree-1 nodes are the two **Gauss points in the interior of the
cell**, not the vertices:

```
the extruded mesh's OWN coordinate element : TensorProductElement(DG1(equispaced), CG1)
what "DG", 1 gives you                     : TensorProductElement(DG1,            DG1)

node z-coords, own coord element : min 0.0000   (cell height 500 m)
node z-coords, DG1 element       : min 105.66   <-- inset by 0.2113 of a cell
```

Neighbouring columns therefore share **no node at all**. Each column fits its own
independent straight line to $z_s(x)$ through its two interior Gauss points, the fits
disagree at the shared facet, and the mesh is torn open along every vertical facet.

**Verification** (180 × 70, `a = 1000 m`, `h_m = 1 m`), measuring $[\![z]\!]$ across
vertical facets — the $z$ coordinate is used because $x$ wraps at the periodic seam and
pollutes the measurement:

| coordinate space | RMS tear | worst tear |
| --- | --- | --- |
| `"DG", 1` (default, spectral) | `2.3e-3 m` | **`2.7e-2 m` = 2.7% of the mountain** |
| `"DG", 1, variant='equispaced'` / own coord element | `1.2e-12 m` | `0.0` ✅ |

The flow was seeing a staircase of ~2.7 cm steps instead of a smooth 1 m hill, with the
steps largest exactly over the peak where $z_s''$ is largest — and each step radiates
its own grid-scale response. This is a *geometric* discontinuity, so no amount of solver
convergence removes it.

`utils.high_dim_mesh_hierarchy` gets away with the same code only because its map
$(x,z)\mapsto(x,0,z)$ is linear and so is reproduced exactly at any set of nodes. The
moment a nonlinear orography is composed with it, the node positions matter.

**Fix:** interpolate the coordinates into a space whose nodes are the cell vertices —
either `variant='equispaced'`, or, more robustly, the mesh's own coordinate element:

```python
coord_elt = m.coordinates.function_space().ufl_element().sub_elements[0]
coord_fs = VectorFunctionSpace(m, coord_elt, dim=dim)
```

### 6.3 Still open: no vertical stabilisation for buoyancy transport

Not a bug so much as a missing term, flagged for completeness. `Vb` is
$DG_h \otimes CG_v$ — **continuous in the vertical** — so `jump(q)` vanishes on `dS_h`
and the upwind facet term in `Nonlinear_buoyancy_Irk` contributes *nothing* on
horizontal facets. Vertical transport of $b$ is therefore a pure centred Galerkin
scheme with no stabilisation at all. The reference `theta_tendency` adds exactly the
term that handles this:

```python
h = avg(CellVolume(mesh))/FacetArea(mesh)
eqn += h**2*c_pen*abs(inner(u('+'), n('+'))) \
       * inner(jump(grad(theta)), jump(grad(q)))*(dS_v + dS_h)   # c_pen = 2**(-7/2)
```

If grid-scale structure survives in $b$ once 6.1 and 6.2 are fixed, this is the next
thing to add.

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
