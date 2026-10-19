'''
Apply the Table-5 edits to ../../../Patch-Preconditioner-for-Boussinesq/ex_article.tex from the Set-T/H run data
(via table5_final.collect()). Prints the numbers it uses. Run once; it refuses to run twice (asserts the old text).
'''
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import table5_final as tf          # noqa: E402
import parse_table5 as pt          # noqa: E402

TEX = os.path.join(HERE, '..', '..', '..', 'Patch-Preconditioner-for-Boussinesq', 'ex_article.tex')

out, hyb = tf.collect()
MESH = tf.MESHES
fs = {m: out[('fs',) + m] for m in MESH if ('fs',) + m in out}
mo = {m: out[('monofix',) + m] for m in MESH if ('monofix',) + m in out}
assert len(fs) == len(MESH) and len(mo) == len(MESH), 'Set T incomplete'
dits = [mo[m]['its'][0] - fs[m]['its'][0] for m in MESH]
per = [(mo[m]['t_med'] / mo[m]['its'][0]) / (fs[m]['t_med'] / fs[m]['its'][0]) for m in MESH]
tot = [mo[m]['t_med'] / fs[m]['t_med'] for m in MESH]
print('its difference mono - fs:', dits)
print('per-iteration cost ratio mono/fs:', [round(x, 3) for x in per])
print('total PCApply ratio mono/fs:', [round(x, 3) for x in tot])
R = {r['tag']: r for r in pt.parse_runs()}
h = {m: R.get(f'H_mono_right_np16_nx{m[0]}_nz{m[1]}') for m in [(20, 20), (40, 40), (80, 80)]}
print('old configuration (FGMRES right):', {m: (r['its_all'], r['reason']) if r else None for m, r in h.items()})

tab = open(os.path.join(HERE, 'table5_tabular.tex')).read().rstrip('\n')
nmin, nmax = tf.ndofs(*MESH[0]), tf.ndofs(*MESH[-1])
pct = lambda x: int(round(100 * (x - 1)))                       # noqa: E731
pmin, pmax = min(pct(x) for x in per), max(pct(x) for x in per)
tmin, tmax = min(pct(x) for x in tot), max(pct(x) for x in tot)
same = all(d == 0 for d in dits)
fewer = [m for m, d in zip(MESH, dits) if d < 0]
its_sentence = ('The two converge in the same number of iterations at every resolution.' if same else
                'The two converge in the same number of iterations, the monolithic cycle needing one fewer at '
                + ', '.join(f'$({a}, {b})$' for a, b in fewer) + '.')
h_its = [h[m]['its'] if h[m] and h[m]['reason'] == 'CONVERGED_RTOL' else None for m in [(20, 20), (40, 40), (80, 80)]]

s = open(TEX).read()
E = []
# 1. timing section: processes of the two tables
E.append((r'''relative to a monolithic alternative. Unless stated otherwise, the
timings below are obtained on $16$ MPI processes on a workstation
with $64\,\mathrm{GB}$ of memory, using the two-dimensional extruded
vertical-slice mesh of Section~\ref{sec:discretisation}.''',
r'''relative to a monolithic alternative. The timings below are obtained
on a workstation with $64\,\mathrm{GB}$ of memory, on $2$ MPI processes
for Table~\ref{tab:rcm} and on $16$ for Table~\ref{tab:monolithic},
using the two-dimensional extruded vertical-slice mesh of
Section~\ref{sec:discretisation}.'''))
# 2. Table 4 caption: the aspect ratio of the runs
E.append((r'''    factorisations. Aspect ratio $\varepsilon \sim 10^{-3}$,
    $\Delta t = 100\,\mathrm{s}$, $\delta = 10^{-4}$.}
    \label{tab:rcm}''',
r'''    factorisations. Aspect ratio $\varepsilon = 1.3\times10^{-2}$
    ($L_x = 300\,\mathrm{km}$, $L_z = 4\,\mathrm{km}$),
    $\Delta t = 100\,\mathrm{s}$, $\delta = 10^{-4}$.}
    \label{tab:rcm}'''))
# 3. the Table-5 paragraph
i0 = s.index(r'''Table~\ref{tab:monolithic} reports, for the same parameter values as''')
i1 = s.index(r'''per level, is left for future tuning.''', i0) + len(r'''per level, is left for future tuning.''')
old_par = s[i0:i1]
hyb_sentence = ''
if h_its[0] and h_its[1]:
    hyb_sentence = (r''' A monolithic cycle whose finest-level smoother acts on the unshifted
operator instead, with $A_\delta$ used only for the patch
factorisations and the coarser levels (the default when PETSc's
multigrid is handed the outer operator, as Firedrake's auxiliary-operator
preconditioner does), is a different
preconditioner: it approximates the inverse of the unshifted operator
rather than $A_\delta^{-1}$, and it needs fewer iterations on the
coarser meshes (''' + f'${h_its[0]}$ and ${h_its[1]}$ at $(20, 20)$ and $(40, 40)$' + r'''), fewer even than an
exact solve with $A_\delta$, but ''' + (r'''it does not converge within $100$
iterations at $(80, 80)$, so it is not robust.''' if not h_its[2] else
f'''its count grows to ${h_its[2]}$ at $(80, 80)$.'''))
new_par = (r'''Table~\ref{tab:monolithic} reports, for the geometry, time step and
shift of Table~\ref{tab:rcm} and five mesh resolutions from
''' + f'${nmin / 1e3:.0f}\\times10^{{3}}$ to ${nmax / 1e6:.1f}\\times10^{{6}}$' + r''' unknowns, the number of
degrees of freedom per patch on the finest level, the cumulative
wall-clock cost of the preconditioner applications (PCApply) over one
FGMRES solve, and the iteration count to a relative residual of
$10^{-8}$, for both approaches with one multigrid cycle per application
and the same smoother. ''' + its_sentence + r'''  This is to be expected.
The pressure space is discontinuous, so its mass matrix is block
diagonal, and every pressure unknown that the divergence couples to the
velocity unknowns of a column-star patch lies in the two columns of
that patch.  Eliminating the pressure from a monolithic patch problem
therefore gives exactly the fieldsplit patch problem on $S$; the
coarse-grid corrections agree for the nested spaces used here; and one
monolithic cycle on $A_\delta$ returns the same velocity as one
fieldsplit cycle.  The two differ only in the pressure, which the Schur
factorisation recovers exactly from the shifted continuity equation and
the monolithic cycle returns with a small residual.  With equal
iteration counts the cost is decided by the patches, which in the
monolithic case also carry the pressure unknowns of their two columns,
$38N_z + 8$ against $30N_z + 8$: an application of the monolithic cycle
is ''' + f'${pmin}$--${pmax}\\%$' + r''' more expensive, and the cumulative PCApply time
of the fieldsplit Schur approach is ''' + (f'lower at every resolution, by ${tmin}$--${tmax}\\%$' if tmin > 0 else
f'within ${max(abs(tmin), abs(tmax))}\\%$ of the monolithic one') + r'''.''' + hyb_sentence)
E.append((old_par, new_par))
# 4. Table 5 caption and tabular
j0 = s.index(r'''    \caption{Fieldsplit Schur-complement preconditioning versus a
    monolithic patch solver applied directly to $A_\delta$.''')
j1 = s.index(r'''    \end{tabular}''', j0) + len(r'''    \end{tabular}''')
new_tab = (r'''    \caption{Fieldsplit Schur-complement preconditioning versus a
    monolithic patch smoother applied to $A_\delta$, with one multigrid
    cycle per application in both (Chebyshev(6) smoothing with the
    column-star patch and RCM-reordered patch factorisations on all
    three levels).  FGMRES with right preconditioning to a relative
    residual of $10^{-8}$; aspect ratio $\varepsilon = 1.3\times10^{-2}$,
    $\Delta t = 100\,\mathrm{s}$, $\delta = 10^{-4}$, $16$ MPI processes.
    $N_x$ is the number of columns of the coarsest multigrid level (the
    finest has $4N_x$); patch DOFs are those of the finest level; times
    are medians of three runs.}
    \label{tab:monolithic}
''' + tab)
E.append((s[j0:j1], new_tab))
# 5. conclusions
E.append((r'''And
the comparison against a monolithic column-patch smoother applied
directly to $A_\delta$ is, in the present untuned configuration,
not in favour of the Schur reduction.
The fieldsplit patches are
about $20\%$ smaller and cheaper per iteration, but the monolithic
patches resolve more of the pressure--velocity coupling per sweep and
converge in fewer FGMRES iterations.''',
r'''And a
monolithic column-patch smoother applied directly to $A_\delta$
converges in the same number of iterations as the Schur reduction, as
it must, since eliminating the pressure from its patches gives the
fieldsplit patches; the fieldsplit patches are about $20\%$ smaller,
so the Schur reduction is the cheaper of the two.'''))
for a, b in E:
    assert s.count(a) == 1, a[:90]
    s = s.replace(a, b)
open(TEX, 'w').write(s)
print(len(E), 'edits applied to', os.path.abspath(TEX))
print(new_par)
