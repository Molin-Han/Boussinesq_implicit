# Lock exchange: Boussinesq vertical slice test case

Notes accompanying `Nonlinear_Boussinesq_Lock_Exchange.py`.

Reference material:

- H. R. Hiester, M. D. Piggott and P. A. Allison, *The impact of mesh adaptivity on the
  gravity current front speed in a two-dimensional lock-exchange*, Ocean Modelling
  **38** (2011) 1–21,
  [doi:10.1016/j.ocemod.2011.01.003](https://doi.org/10.1016/j.ocemod.2011.01.003).
  This is the test case requested.
- C. Härtel, E. Meiburg and F. Necker, *Analysis and direct numerical simulation of the
  flow at a gravity-current head. Part 1*, J. Fluid Mech. **418** (2000) 189–212. The
  configuration Hiester et al. adopt, and the source of the reference front speeds.
- L. Yue, O. B. Fringer et al., *An unstructured grid, nonhydrostatic, generalized
  vertical coordinate ocean model*, [arXiv:2109.07467](https://arxiv.org/abs/2109.07467),
  section 5.3. Restates the same configuration in full, and is open access — see
  section 2 on why this matters.

---

## 1. The test case

A closed, two-dimensional channel of length $L = 0.8$ m and height $D = 0.1$ m, with a
lock gate at $x_c = L/2$. Dense fluid fills the left half, light fluid the right half,
and everything is at rest. At $t = 0$ the gate is removed; the dense fluid runs
rightwards along the floor and the light fluid leftwards along the ceiling, and the two
counter-propagating gravity currents shear against each other, rolling up into
Kelvin–Helmholtz billows along the interface.

The non-dimensional density difference is $\Delta\rho/\rho_0 = 10^{-3}$, giving a reduced
gravity

$$g' = g\,\frac{\Delta\rho}{\rho_0} = 0.01\ \mathrm{m\,s^{-2}},$$

a buoyancy velocity and a gravity-current time scale built on the **half** depth,

$$u_b = \sqrt{g' D/2} = 0.02236\ \mathrm{m\,s^{-1}}, \qquad
T_b = \sqrt{D/(2g')} = 2.236\ \mathrm{s},$$

and, with molecular viscosity $\nu = 10^{-6}\ \mathrm{m^2 s^{-1}}$ and **no** scalar
diffusivity ($\kappa = 0$, so $Sc = \infty$), a Grashof number

$$Gr = \left(\frac{u_b\,(D/2)}{\nu}\right)^{2} = 1.25\times10^{6}.$$

The boundary conditions are the distinctive part of this configuration: no normal flow
everywhere, but **no slip on the bottom** and **free slip (rigid lid) on the top**. A
single run therefore delivers both a no-slip and a free-slip gravity current, which is
the whole point — the two have measurably different front speeds and stress the
discretisation in different ways.

**The diagnostic** is the non-dimensional front speed, i.e. the Froude number
$Fr = u_g/u_b$ of each current, obtained by linear regression of the front position
against time. Both fronts accelerate out of the lock, then travel at a roughly steady
speed, and hit the end walls at $0.4$ m from the lock; the fit is taken over the window
in which a front has travelled between $0.2$ m and $0.3$ m.

Reference values, from the Härtel et al. (2000) DNS:

| current | wall | $Fr$ |
| --- | --- | --- |
| light, leftward, along the ceiling | free slip | **0.675** |
| dense, rightward, along the floor | no slip | **0.574** |

---

## 2. Provenance of the numbers — please check this section

> **I could not read the Hiester et al. paper directly.** It is paywalled and every
> route I tried (the DOI, the ScienceDirect page, the Daneshyari preview PDF) returned
> a 403, a 504 or a redirect I could not follow. Everything in section 1 is taken from
> **Yue et al. (2021), section 5.3**, which reproduces this configuration, cites Hiester
> et al. (2011) as one of the codes that runs it, and states each parameter explicitly;
> the domain size and the Grashof number were independently corroborated by secondary
> sources. You have access to the paper and I do not, so this section is the one worth
> proofreading hardest.

Confirmed by an explicit statement in a source I read:

| quantity | value | source |
| --- | --- | --- |
| channel length $L$ | 0.8 m | Yue et al. §5.3, and secondary sources |
| channel height $D$ | 0.1 m | Yue et al. §5.3 |
| $\Delta\rho/\rho_0$ | $10^{-3}$ | Yue et al. §5.3 |
| $g'$ | 0.01 m s⁻² | Yue et al. §5.3 |
| $u_b = \sqrt{g'D/2}$ | 0.022 m s⁻¹ | Yue et al. §5.3 |
| $T_b = \sqrt{D/2g'}$ | 2.24 s | Yue et al. §5.3 |
| $\nu$ | $10^{-6}$ m² s⁻¹ | Yue et al. §5.3, attributed to Hiester et al. |
| $\kappa$ | 0 | Yue et al. §5.3, attributed to Hiester et al. |
| $Gr$ | $1.25\times10^{6}$ | attributed to Hiester et al. |
| no slip bottom / free slip lid | — | Yue et al. §5.3, Fig. 6 |
| $Fr$ = 0.675 free slip, 0.574 no slip | — | Yue et al. Table 3, quoting Härtel et al. |
| fit window: 0.2 m to 0.3 m travelled | — | secondary source describing Hiester's method |
| fronts reach the end walls at 0.4 m | — | same |

**Open points I could not confirm, and what the script assumes:**

1. **End-wall conditions.** The script uses free slip on the two vertical end walls.
   They are only felt once a front arrives, which is after the measuring window closes,
   so this should not affect the diagnostic — but it is an assumption.
2. **Sharpness of the initial interface.** The script defaults to a **discontinuous**
   lock gate, which is what the description implies. `--interface_width` applies a
   $\tanh$ smoothing if the paper in fact uses one. Worth checking: a discontinuous
   initial condition in a DG space without a limiter will produce some over- and
   undershoot near the gate.
3. **Which side is dense.** Immaterial by symmetry, but the script puts dense fluid on
   the **left**, so that the no-slip current is the rightward one and the free-slip
   current the leftward one — matching the orientation of Yue et al.'s figure and hence
   their Table 3 labelling.
4. **End time and time step.** Not stated in what I could read. The script's defaults
   are mine (section 4.7), chosen so the fit window is comfortably covered.
5. **Resolution.** Hiester et al.'s point is mesh adaptivity, so there is no single
   "paper resolution" to reproduce. The default $256\times64$ is the resolution of
   Yue et al.'s fixed-grid case Z64, which is a reasonable fixed-mesh reference point.

---

## 3. The mathematical problem in the Boussinesq setting

With kinematic pressure $p$ and buoyancy $b$, and **no background stratification**
($N = 0$), non-rotating:

$$\frac{\partial u}{\partial t} + (u\cdot\nabla)u = -\nabla p + b\,\hat k + \nu\nabla^2 u,$$

$$\frac{\partial b}{\partial t} + (u\cdot\nabla)b = 0,$$

$$\nabla\cdot u = 0,$$

on $0 \le x \le L$, $0 \le z \le D$, with $u\cdot n = 0$ on all four walls, $u = 0$ on
the floor and zero tangential stress on the ceiling and end walls.

**Initial condition:** $u = 0$, $p = 0$, and

$$b(x, z, 0) = \begin{cases} -g'/2, & x < x_c \quad \text{(dense)},\\[2pt]
+g'/2, & x > x_c \quad \text{(light)}. \end{cases}$$

Note the sign convention: `utils.Nonlinear_velocity_Irk` carries the buoyancy as
$+b\,\hat k$, so negative $b$ is dense fluid and sinks. The buoyancy jump across the gate
is exactly $g'$.

Because $N = 0$, the linear stratification term $N^2(u\cdot\hat k)$ in
`utils.Nonlinear_buoyancy_Irk` is switched off by passing `N2=Constant(0.0)`, and the
buoyancy equation is pure advection. All of the buoyancy lives in the prognostic field.

**Contrast with the mountain case.** There, the initial state is an exact steady
solution and the orography is the sole, tiny forcing — the whole problem is a $10^{-3}$
perturbation. Here the initial state is maximally *out* of balance: it is a finite step
in $b$ with nothing holding it up, and the entire flow is the collapse of that step. The
nonlinear advection term carries the solution from the first step onwards, and the
available potential energy released,

$$E_p(0) = -\!\int b\,z\ \mathrm{d}x = 0 \quad\text{by symmetry},$$

converts into kinetic energy as the currents run. That makes it a much stiffer test of
the advection discretisation than the mountain, and a useful complement to it.

---

## 4. What the script implements

`Nonlinear_Boussinesq_Lock_Exchange.py` keeps the skeleton of
`Nonlinear_Boussinesq_Mountain.py` unchanged: `utils` residual forms, the 3D-embedded
slice mesh, `MeshHierarchy` + `ExtrudedMeshHierarchy`, Irksome `GaussLegendre(1)`
(implicit midpoint), and the `HDivSchurPC(IRKAuxiliaryOperatorPC)` shift preconditioner
with the $p \to p - \delta^{-1}\nabla\cdot u$ elimination, fieldsplit-Schur, a mass solve
on the pressure block and `mg` / `ASMStarPC` on the velocity–buoyancy block. Element
degrees and solver parameters are untouched.

The deviations, in rough order of how much thought they needed:

### 4.1 A viscous term — symmetric interior penalty

The mountain test case is inviscid; here $\nu$ *defines* the Grashof number, and the
no-slip front speed is set by the bottom boundary layer, so viscosity cannot be dropped.

The velocity lives in $\mathrm{RT}\times DG$ and is only $H(\mathrm{div})$-conforming:
the normal component is continuous across facets but the tangential one jumps. A viscous
term therefore needs a DG-style treatment of those tangential jumps. `viscous_form`
assembles the symmetric interior penalty (SIPG) form of $-\nu\nabla^2 u$,

$$\nu\!\int\!\nabla w : \nabla u
\;-\;\nu\!\int_\Gamma\!\Big(\{\nabla w\} : [\![u\otimes n]\!] + [\![w\otimes n]\!] : \{\nabla u\}\Big)
\;+\;\nu\!\int_\Gamma \sigma\,[\![w]\!]\cdot[\![u]\!],$$

summed over `dS_v` and `dS_h` separately.

**No slip has to be imposed weakly.** This is the subtle point. A
`DirichletBC(W.sub(0), 0, "bottom")` on an extruded RT space constrains the *facet*
degrees of freedom of the bottom, which carry the **normal** (vertical) velocity only —
the wall-tangential velocity is a cell-interior degree of freedom and the strong
boundary condition never touches it. The strong BC therefore gives $u\cdot n = 0$ and
nothing more. The tangential no-slip condition is added by a Nitsche term on `ds_b`,

$$-\nu\!\int_{\Gamma_b}\!\Big(\nabla w : (u\otimes n) + (w\otimes n) : \nabla u\Big)
\;+\;\nu\!\int_{\Gamma_b}\!\sigma\, w\cdot u .$$

Boundaries with **no** such term get the natural condition, zero tangential stress,
which is exactly free slip. So the lid and the end walls need no code at all — free slip
is what you get by saying nothing, and that is why `no_slip_measures = (ds_b,)` is the
only entry.

**The penalty uses `CellDiameter`, not a precomputed constant.** $\sigma = C p^2/h$ with
$h$ the UFL cell diameter, so that when the auxiliary operator of the preconditioner is
coarsened onto the multigrid hierarchy the penalty rescales with the coarse cells. A
frozen fine-grid $\sigma$ would be far too small on the coarse levels and lose
coercivity there. $C$ is exposed as `--ip_penalty`, default 10.

Cell aspect ratio at the default resolution is only $\Delta x/\Delta z = 2$, so using the
isotropic cell diameter on both `dS_v` and `dS_h` over-penalises the vertical facets by
a factor of about 2.2. With $\nu = 10^{-6}$ the viscous block is a small perturbation of
the operator, so this costs essentially nothing; if the aspect ratio is ever pushed hard,
direction-specific penalties would be the thing to add.

### 4.2 Closed channel, no orography

`PeriodicIntervalMesh` becomes `IntervalMesh` — the lock exchange has end walls. The two
ends of the interval keep the markers `1` ($x = 0$) and `2` ($x = L$) through the
extrusion and the 3D embedding, and become the side walls of the slice, so the boundary
conditions are

```python
bcs = [DirichletBC(W.sub(0), 0, "top"),  DirichletBC(W.sub(0), 0, "bottom"),
       DirichletBC(W.sub(0), 0, 1),      DirichletBC(W.sub(0), 0, 2)]
```

There is no terrain, so `mountain_mesh_hierarchy` is not needed and the plain
`utils.high_dim_mesh_hierarchy` embedding $(x,z)\mapsto(x,0,z)$ is used. That map is
linear and hence exact at any set of nodes, which is precisely why it can use the default
DG variant where the mountain script has to ask for `variant='equispaced'` (see §6.2 of
the mountain notes).

### 4.3 The initial condition is trivial — no Leray projection

The mountain script needs a discrete Leray projection of its uniform wind, because
Irksome's `bc_type="DAE"` constrains only the stage values: with the implicit midpoint
rule that means $\operatorname{div} u^{n+1} = -\operatorname{div} u^{n}$ and
$u^{n+1}\!\cdot n|_\partial = -u^{n}\!\cdot n|_\partial$, so an initial state that is not
solenoidal or not compatible with the boundary conditions is never cleaned up — it just
flips sign every step.

Here the fluid starts **at rest**. $u_0 = 0$ satisfies both conditions exactly, and the
projection machinery is simply not required. Confirmed at run time:
$\lVert\operatorname{div}u_0\rVert_{L^2} = 0$ exactly.

### 4.4 Front position diagnostic

The natural definition — chase a contour of $b$ — becomes ill-defined once the billows
start shedding isolated patches, and needs point evaluation, which is awkward in
parallel. Instead:

> In a thin strip of thickness $h_p$ along a wall, the current is a single tongue
> attached to the end wall it came from. So the **area** of dense fluid in the bottom
> strip, divided by the strip's area and multiplied by $L$, is the $x$-coordinate of the
> nose of the bottom current; and symmetrically from the right for the top one.

This is two `assemble` calls of an indicator function — cheap, collective, no point
location — so the fronts are tracked **every time step**, not every dump. That matters:
the Froude number is a regression over a window covering only about a tenth of the run,
and it wants every sample it can get.

Two details are needed to make it exact rather than merely close:

- the strip is **snapped to a whole number of layers**, so that it ends on a facet. An
  indicator that cuts through the middle of a cell is only integrated to quadrature
  accuracy;
- the normalisation is by the **measured** strip area, not by $L\,h_p$, which removes
  whatever bias is left.

Without those, the front position at $t=0$ came out as $0.4130$ m instead of $0.4$ m — a
systematic offset of a quarter of a cell that would have contaminated the fit. With
them, $t=0$ gives $0.4000000000000025$ m.

The strip thickness is `--front_depth_frac`, default $0.1D$ (6 layers, $h_p = 9.4$ mm at
the default $n_z = 64$). This is the one diagnostic parameter I would sanity-check
against the paper's own definition if it states one: too thin and the strip misses the
raised nose of the no-slip head, too thick and it starts measuring the body of the
current rather than its front.

The Froude numbers are fitted at the end of the run by `froude_number`, printed against
the DNS values with a percentage error, and the whole time series is written to
`Nonlinear_Boussinesq_lock_exchange_fronts.txt` for replotting.

### 4.5 Energetics

`KE` $= \tfrac12\int|u|^2$ and `PE` $= -\int b\,z$ are printed at every dump. Since
$E_p(0) = 0$ by symmetry, `PE` *is* the released potential energy, and with $\kappa = 0$
the gap between $|{\rm PE}|$ and `KE` is dissipation — viscous plus numerical. It is a
cheap and quite sensitive health check on the advection scheme.

### 4.6 Shift parameter

$\delta$ has units of s/m², so as in the mountain case the tuned unit-box value
$10^{-4}$ does not transfer. There is no mean flow here, so the buoyancy velocity plays
the role of $U$:

$$\delta = \frac{10^{-4}}{u_b L} \approx 5.6\times10^{-3}\ \mathrm{s\,m^{-2}},$$

exposed as `--shift_scale` (the dimensionless group $\delta u_b L$) with a raw override
via `--shift`. This has **not** been tuned at these scales — it is transferred by
dimensional analysis only, exactly as flagged at the end of the mountain notes. It is
the first parameter to scan if the KSP misbehaves.

### 4.7 Defaults

| quantity | default | reasoning |
| --- | --- | --- |
| columns × layers | 256 × 64 | `--nx 64 --refinement 2` (64·2² = 256), `--nz 64`; the resolution of Yue et al.'s fixed-grid case Z64 |
| $\Delta x$, $\Delta z$ | 3.125 mm, 1.5625 mm | |
| `--dt` | 0.02 s | advective Courant $u_b\Delta t/\Delta x \approx 0.14$, close to the 0.18 of the explicit reference codes. Implicit midpoint is unconditionally stable, so this is an accuracy choice, not a stability one |
| `--tmax` | 25 s $\approx 11 T_b$ | the fit window closes at about $t = 23$ s; the fronts reach the end walls at about 26 s (top) and 31 s (bottom) |
| `--dumpt` | 0.5 s | 50 VTK dumps; the mountain script dumps every step, which would be 1250 here |
| `--atol` | $10^{-30}$ | as in the mountain script: $g' = 10^{-2}$ and velocities are $O(10^{-2})$, so absolute residuals are tiny and the KSP must be driven by `rtol` |

---

## 5. Verification so far

### 5.1 What was checked

Two runs, both completed cleanly (exit 0). Output files were deleted afterwards.

**(a) Tiny smoke test** — $16\times8$, `--refinement 1`, `dt = 0.05`, 2 steps:

| check | result |
| --- | --- |
| $\lVert\operatorname{div}u_0\rVert_{L^2}$ | `0.0` exactly |
| initial buoyancy range | `[-0.005, 0.005]` = $\pm g'/2$ ✅ |
| Newton, step 1 | `3.68e-4 → 4.08e-8 → 2.54e-14` (2 iterations) |
| linear solves | converge by `rtol` in 4–5 iterations |
| KSP iterations per step | 6.5 |
| `max|uy|` | `0.0` |

**(b) Coarse physics run** — $64\times16$, `dt = 0.1`, to $t = 12$ s $= 5.4\,T_b$:

| check | result |
| --- | --- |
| front position at $t=0$ | `0.4000000000000025` m, both fronts ✅ |
| direction of travel | dense front rightwards along the floor, light front leftwards along the ceiling ✅ |
| KE gain vs. PE release at $t=12$ s | `3.42e-6` vs. `3.72e-6`, i.e. 92% converted, 8% dissipated |
| `max|uy|` throughout | `0.0` |
| max Courant | 0.22 |
| linear solves | 4–6 iterations |

`max|uy| = 0.0` is the health check worth keeping: nothing in this test case forces the
$y$ direction, so any growth there would mean the embedding of the slice in 3D had gone
wrong.

Fitting that coarse run over a *narrowed* window (0.05–0.12 m travelled, since 12 s is
not long enough to reach 0.2 m) exercises the whole regression path and gives:

| current | $Fr$ | DNS | error |
| --- | --- | --- | --- |
| no slip (bottom) | 0.484 | 0.574 | −15.7 % |
| free slip (top) | 0.736 | 0.675 | +9.1 % |

Both errors are the expected sign and size for $64\times16$. $\Delta z = 6.25$ mm leaves
the bottom boundary layer completely unresolved, which is exactly the mechanism that
limits the no-slip head speed — Yue et al. make the same observation about their own
coarse cases. The free-slip front reads high because at 0.05–0.12 m the current is still
in its acceleration phase, and because the wall strip is $0.125D$ thick on that mesh.

### 5.2 What has *not* been checked

> **The default $256\times64$ configuration has not been run.** I have no evidence about
> how close the Froude numbers land to 0.574 / 0.675 at the intended resolution, and the
> section 5.1 numbers should not be read as a prediction. Convergence of $Fr$ under
> refinement is the actual acceptance test for this script, and it is still to do.

Also outstanding:

- **No visual check of the Kelvin–Helmholtz billows.** The reference papers make a point
  of the interface roll-up being captured; nobody has looked at the `.pvd` yet.
- **No check of the initial discontinuity's over/undershoot.** A sharp gate in a DG space
  with no limiter will ring. Whether it matters for the front speed is unknown.
- **No parallel run.** All diagnostics are written to be collective, and the only
  rank-local reductions (`global_max`, the initial buoyancy range) go through
  `allreduce`, but this is untested.
- **`--shift_scale` untuned**, see §4.6.
- The mountain notes' §6.3 caveat applies here unchanged: `Vb` is continuous in the
  vertical, so `jump(q)` vanishes on `dS_h` and vertical transport of $b$ is an
  unstabilised centred scheme. This test case advects a *discontinuous* $b$ field
  vigorously in the vertical, so if that missing stabilisation term is going to show up
  anywhere, it will show up here.

---

## Running it

```bash
# defaults: 256 x 64, dt = 0.02 s, to t = 25 s
python Nonlinear_Boussinesq_Lock_Exchange.py

# a quick look at the physics
python Nonlinear_Boussinesq_Lock_Exchange.py --nx 16 --nz 16 --dt 0.1 --tmax 12 --dumpt 1

# inviscid variant: drops the interior penalty term and makes the floor free slip too
python Nonlinear_Boussinesq_Lock_Exchange.py --nu 0.0

# see all options
python Nonlinear_Boussinesq_Lock_Exchange.py --help
```

The run prints the fitted Froude numbers and their error against the DNS at the end, and
writes the front time series to `Nonlinear_Boussinesq_lock_exchange_fronts.txt`. The
buoyancy field in `Nonlinear_Boussinesq_lock_exchange.pvd` should show the two currents
and the billow train along the interface.
