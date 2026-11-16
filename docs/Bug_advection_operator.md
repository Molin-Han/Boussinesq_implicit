# Bug 1: the nonlinear advection operator in `utils.py`

What was wrong, why it is wrong, and why it showed up as cell-to-cell jumps in the
velocity field of the mountain test case.

---

## 1. The one-line summary

```python
eqn -= inner(div(outer(u, w)), u) * dx     # WAS  -- wrong operator
eqn -= inner(div(outer(w, u)), u) * dx     # NOW  -- correct
```

The two arguments of `outer` were transposed. This is not a sign error or a factor of
two: it assembles a **different differential operator**, one which is (cellwise) a pure
gradient plus a facet-supported remainder. The physical advection term was effectively
missing, and what replaced it was grid-scale noise.

---

## 2. The UFL conventions that make this a real difference

Two conventions matter, and together they mean `outer(u,w)` and `outer(w,u)` are *not*
interchangeable under `div`:

1. **`outer`**: $\texttt{outer}(a,b)_{ij} = a_i\,\overline{b_j}$ — the first argument
   carries the *row* index.
2. **`div` of a rank-2 tensor**: UFL contracts over the **last** index,
   $\texttt{div}(A)_i = \partial_j A_{ij}$.

So

$$\texttt{div(outer(w,u))}_i = \partial_j\!\left(w_i u_j\right), \qquad
\texttt{div(outer(u,w))}_i = \partial_j\!\left(u_i w_j\right).$$

These are the divergences of a tensor and of its transpose. Nothing forces them to
agree.

---

## 3. Deriving the correct weak form

The term we want in the momentum equation is $(u\cdot\nabla)u$, tested against $w$:

$$N(u;w) \;=\; \int_\Omega w\cdot(u\cdot\nabla)u \,\mathrm{d}x
\;=\; \int_\Omega w_i\,u_j\,\partial_j u_i \,\mathrm{d}x .$$

The velocity lives in an H(div) space, so it is only **normal**-continuous: the
tangential components jump across facets and $u_j\partial_j u_i$ is not integrable in
the naive sense. The standard remedy is to integrate by parts cell by cell and then
replace the resulting inter-element flux by an upwind numerical flux.

**Cellwise integration by parts.** On a cell $K$,

$$\int_K w_i u_j \partial_j u_i \,\mathrm{d}x
= -\int_K \partial_j\!\left(w_i u_j\right) u_i \,\mathrm{d}x
+ \oint_{\partial K} (w\cdot u)\,(u\cdot n)\,\mathrm{d}s .$$

The volume term is precisely $-\int_K \texttt{div(outer(w,u))}\cdot u$, by convention (2)
above. **This is where `outer(w, u)` comes from.**

**Summing and upwinding.** Adding over all cells and replacing the interior-facet
contributions by the upwind flux $\tilde u_n = \tfrac12\!\left(u\cdot n + |u\cdot n|\right)$
gives the discrete operator

$$N_h(u;w) = -\int_\Omega \texttt{div(outer(w,u))}\cdot u\,\mathrm{d}x
\;+\; \sum_F \int_F [\![w]\!]\cdot\!\left(\tilde u_n^+ u^+ - \tilde u_n^- u^-\right)\mathrm{d}S,$$

which is exactly what the code assembles once the fix is in:

```python
eqn -= inner(div(outer(w, u)), u) * dx
eqn += dot(jump(w), unn('+')*u('+') - unn('-')*u('-')) * (dS_v + dS_h)
```

This matches `slice_utils.u_tendency` in Cotter's repository (its
`vector_invariant=False` branch), which was the cross-check that first flagged the
discrepancy.

**A useful simplification.** Expanding the correct volume integrand,

$$\texttt{div(outer(w,u))}\cdot u = u_i\,\partial_j(w_i u_j)
= u_i\,u_j\partial_j w_i + (u\cdot w)\,\nabla\!\cdot\! u
= u\cdot(u\cdot\nabla)w + (u\cdot w)\,\nabla\!\cdot\! u .$$

Here the second term vanishes **pointwise**, not just weakly: $\nabla\cdot$ maps
$RT_k$ onto $DG_{k-1}$, which is exactly the pressure space, so enforcing
$\int \phi \,\nabla\!\cdot\! u = 0$ for all $\phi \in DG_{k-1}$ forces
$\nabla\!\cdot\! u \equiv 0$. (Measured in the runs: $\|\nabla\!\cdot\! u\| \sim 2\times10^{-9}$,
i.e. roundoff at these scales.) So the correct volume term is just
$-\int u\cdot(u\cdot\nabla)w$, the transpose-derivative form you would write by hand.

---

## 4. What the buggy form actually computed

Now expand the version that was in the code. Write $K = \tfrac12|u|^2$ for the kinetic
energy density, so $|u|^2 = 2K$:

$$\texttt{div(outer(u,w))}\cdot u
= u_i\,\partial_j\!\left(u_i w_j\right)
= u_i\,w_j\partial_j u_i + u_i u_i\,\partial_j w_j
= (w\cdot\nabla)K + 2K\,\nabla\!\cdot\! w .$$

So the buggy volume term contributed
$-\int_\Omega \left[(w\cdot\nabla)K + 2K\,\nabla\!\cdot\! w\right]\mathrm{d}x$.
Integrating the first piece by parts on each cell,

$$-\int_K (w\cdot\nabla)K \,\mathrm{d}x
= \int_K K\,\nabla\!\cdot\! w \,\mathrm{d}x - \oint_{\partial K} K\,(w\cdot n)\,\mathrm{d}s,$$

and therefore, cell by cell,

$$-\int_K \texttt{div(outer(u,w))}\cdot u \,\mathrm{d}x
= \underbrace{-\int_K K\,\nabla\!\cdot\! w \,\mathrm{d}x}_{\text{(a)}}
\;\underbrace{-\;\oint_{\partial K} K\,(w\cdot n)\,\mathrm{d}s}_{\text{(b)}} .$$

Summing over cells and using that $w\cdot n$ **is** continuous across facets (H(div))
while $K$ is **not** (tangential jumps), term (b) collapses onto the facets:

$$\boxed{\;-\int_\Omega \texttt{div(outer(u,w))}\cdot u \,\mathrm{d}x
= -\int_\Omega K\,\nabla\!\cdot\! w\,\mathrm{d}x
\;-\; \sum_F \int_F [\![K]\!]\,(w\cdot n)\,\mathrm{d}S \;}$$

Read the two pieces:

**(a) is a pure gradient.** Compare it with the pressure term in the momentum equation,
which is `- div(w)*p*dx`. Term (a) is *identical in structure* with $p$ replaced by $K$.
In an incompressible flow a pure gradient is absorbed entirely into the pressure — it
exerts no dynamical influence whatsoever. As far as the physics is concerned, the
advection volume term had been **deleted**, and the solver silently compensated by
shifting $p \mapsto p - K$.

**(b) is grid-scale garbage.** What survives is a term supported *only on facets*,
proportional to the jump in kinetic energy. It has no continuum counterpart; it is a
forcing that lives entirely at cell interfaces.

**And the facet terms are now unbalanced.** The upwind flux terms in the equation were
never touched — they are still the correct numerical fluxes for the *true* advection
operator, derived to balance the volume term from §3. Pairing them with a volume term
that is a pure gradient leaves a residual that is, once again, concentrated on facets.

So the mechanism behind the symptom is:

> the physically meaningful part of the operator was absorbed into the pressure, and
> everything that remained of it lived on cell interfaces.

That is why the corruption appeared as **jumps across cells** in the velocity, rather
than as a smooth error — and why the solver converged perfectly while doing it. It was
solving the wrong equations to tight tolerance, exactly as reported.

---

## 5. Numerical verification

Take a smooth, **continuous**, divergence-free field with $u\cdot n = 0$ on top and
bottom and periodic in $x$: from the streamfunction
$\psi = \sin(2\pi x/L)\sin(\pi z/H)$, set $u = (\partial_z\psi,\,0,\,-\partial_x\psi)$.
Because $u$ and $w$ are continuous, every facet term vanishes and every boundary term
vanishes, so the weak form must reproduce the strong form *exactly*.

| form | rel. difference from $\int w\cdot(u\cdot\nabla)u$ |
| --- | --- |
| `-inner(div(outer(w,u)), u)*dx` — correct | $2.6\times10^{-14}$ ✅ |
| `-inner(div(outer(u,w)), u)*dx` — was in `utils.py` | $7.4\times10^{0}$ ❌ |

Off by 737%, i.e. not the same operator by any margin.

End-to-end in the mountain test case, with a **direct LU solve** so that the
preconditioner is entirely out of the picture (180 columns, 35 layers, $\Delta t = 25$ s):

| advection form | max\|w\| at $t=500$ s | RMS $[\![u]\!]/\lvert u\rvert$ |
| --- | --- | --- |
| `outer(u,w)` (buggy) | $1.30\times10^{-2}$ m/s, drifting | $9.6\times10^{-5}$ |
| `outer(w,u)` (fixed) | $4.88\times10^{-3}$ m/s, steady | $5.1\times10^{-6}$ |

The grid-noise measure drops by a factor of ~19, the solution becomes steady instead of
drifting, and the amplitude lands on linear theory: the analytic surface value is
$U\max|\mathrm{d}z_s/\mathrm{d}x| = 10 \times \tfrac{3\sqrt3}{8}\tfrac{h_m}{a}
= 6.5\times10^{-3}$ m/s. Verified identical through the shift preconditioner
($4.8817\times10^{-3}$) and steady out to $t = 3000$ s.

---

## 6. Why it stayed hidden

The bugged term is quadratic in $u$. Every earlier test in this folder ran with
`--U_mean 0` and a buoyancy perturbation of amplitude $\sim 3\times10^{-4}$, so $|u|$
was tiny and the whole advection term was negligible against the linear terms — those
problems are effectively linear Boussinesq, and a wrong quadratic term contributes
nothing measurable.

The mountain test case is the first one in this folder with a **10 m/s mean flow**. There
the advection term is not a correction, it is the term that transports the entire
mountain wave. The bug went from invisible to dominant in one test case.

---

## 7. Where the fix was applied

All four sites in `utils.py` had the same transposition:

| function | line | note |
| --- | --- | --- |
| `Nonlinear_velocity_Irk` | ~205 | the one used by the mountain script |
| `Nonlinear_velocity` | ~175 | non-Irksome nonlinear version |
| `LB_velocity` | ~118, ~126 | the two `U_mean` mean-flow terms |

Since `utils.py` is shared, this changes `Nonlinear_Boussinesq_Irk_SC.py`,
`Nonlinear_Boussinesq_slice.py` and the `LB_*` mean-flow paths. Per §6 the effect there
should be negligible, but results are no longer byte-identical — worth re-running
anything already published.

---

## 8. The general lesson

`div(outer(a, b))` reads symmetrically but is not. Whenever you write a weak advection
term this way, check it against the strong form on a smooth continuous field, where all
facet and boundary terms drop out and the identity has to hold to machine precision.
It is a five-line test and it would have caught this immediately:

```python
u = Function(VectorFunctionSpace(mesh, "CG", 3, dim=3)).interpolate(u_exact)  # continuous
w = TestFunction(u.function_space())
strong  = assemble(inner(w, dot(grad(u), u))*dx)
weak    = assemble(-inner(div(outer(w, u)), u)*dx)
assert np.allclose(strong.dat.data, weak.dat.data)
```
