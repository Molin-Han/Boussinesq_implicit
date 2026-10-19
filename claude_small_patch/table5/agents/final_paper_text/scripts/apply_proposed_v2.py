"""Apply the referee's proposed replacements to a COPY of ex_article.tex (scratch only).
usage: python -I apply_proposed_v2.py SRC_TEX DST_TEX OUT_JSON"""
import json, sys
src, dst, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
t = open(src).read()
R = []
def add(name, old, new): R.append((name, old, new))

add("R1 RCM paragraph: PCApply is not a single application",
r"""Table~\ref{tab:rcm} reports the patch matrix bandwidth and the
wall-clock cost of a single preconditioner application (PCApply),
with and without RCM reordering, for a fixed aspect ratio""",
r"""Table~\ref{tab:rcm} reports the patch matrix bandwidth and the
wall-clock time spent in the preconditioner (PCApply) during one GMRES
solve, over its $14$ applications, with and without RCM reordering,
for a fixed aspect ratio""")
add("R2 RCM paragraph: header name",
r"""(a least-squares fit of $\mathrm{Time(PCApply)}$ against $N_z$""",
r"""(a least-squares fit of PCApply against $N_z$""")
add("R3 Table 4 caption",
r"""    \caption{Effect of RCM reordering on the column-patch LU
    factorisations. Aspect ratio $\varepsilon = 1.3\times10^{-2}$""",
r"""    \caption{Effect of RCM reordering on the column-patch LU
    factorisations. PCApply is the wall-clock time of the $14$
    preconditioner applications of one GMRES solve on $2$ MPI
    processes, including the numerical refactorisation of the patch
    matrices. Aspect ratio $\varepsilon = 1.3\times10^{-2}$""")
add("R8 unknowns range",
r"""$32\times10^{3}$ to $2.1\times10^{6}$ unknowns""",
r"""$3.2\times10^{4}$ to $2.05\times10^{6}$ unknowns""")
add("R4 Table-5 paragraph: what PCApply contains; 'That they agree'",
r"""degrees of freedom per patch on the finest level, the cumulative
wall-clock cost of the preconditioner applications (PCApply) over one
FGMRES solve, and the iteration count to a relative residual of
$10^{-8}$, for both approaches with one multigrid cycle per application
and the same smoother.  The two converge in the same number of
iterations, except at $(80, 80)$, where the monolithic cycle needs one
fewer.  This is to be expected.""",
r"""degrees of freedom per patch, the wall-clock time spent in the
preconditioner (PCApply) during one FGMRES solve, which includes the
numerical refactorisation of the patch and coarse-grid matrices, and
the iteration count to a relative residual of
$10^{-8}$, for both approaches with one multigrid cycle per application
and the same smoother.  The two converge in the same number of
iterations, except at $(80, 80)$, where the monolithic cycle needs one
fewer.  That they agree is to be expected.""")
add("R5 Table-5 paragraph: coarse grid, Chebyshev interval, cost sentence",
r"""therefore gives exactly the fieldsplit patch problem on $S$; the
coarse-grid corrections agree for the nested spaces used here; and one
monolithic cycle on $A_\delta$ returns the same velocity as one
fieldsplit cycle.  The two differ only in the pressure, which the Schur
factorisation recovers exactly from the shifted continuity equation and
the monolithic cycle returns with a small residual.  The cost is
therefore decided by the patches, which in the monolithic case also
carry the pressure unknowns of their two columns, $38N_z + 8$ against
$30N_z + 8$: an application of the monolithic cycle is $11$--$17\%$
more expensive, and the cumulative PCApply time of the fieldsplit Schur
approach is lower by $11$--$17\%$ at every resolution except
$(80, 80)$, where the monolithic cycle saves an iteration and the
difference is $1\%$.""",
r"""therefore gives exactly the fieldsplit patch problem on $S$; the
coarse-grid corrections agree in the velocity, because the spaces are
nested and the cells affine, so that the rediscretised coarse operators
are the Galerkin ones and the divergence commutes with the
prolongation; and, since PETSc's eigenvalue estimates give both
smoothers the same Chebyshev interval ($[0.2, 2.2]$ to four digits),
one monolithic cycle on $A_\delta$ returns the same velocity as one
fieldsplit cycle.  The two differ only in the pressure, which the
Schur factorisation recovers exactly from the shifted divergence
equation~\eqref{eq:shifted-div} and the monolithic cycle returns with
a small residual.  The difference in cost therefore comes from the
larger systems the monolithic cycle works with: its patches also carry
the pressure unknowns of their two columns, $38N_z + 8$ against
$30N_z + 8$ unknowns, and its smoother, residuals and grid transfers
act on the pressure as well.  An application of the monolithic cycle
is $11$--$17\%$ more expensive, and its cumulative PCApply time is
$11$--$17\%$ higher than that of the fieldsplit Schur approach at every
resolution except $(80, 80)$, where the monolithic cycle saves an
iteration and the difference is $1\%$.""")
add("R6 second paragraph",
r"""A monolithic cycle whose finest-level smoother acts on the unshifted
operator instead, with $A_\delta$ used only for the patch
factorisations and the coarser levels (the default when PETSc's
multigrid is handed the outer operator, as Firedrake's auxiliary-operator
preconditioner does), is a different
preconditioner: it approximates the inverse of the unshifted operator
rather than $A_\delta^{-1}$, and it needs fewer iterations on the
coarser meshes ($6$ and $8$ at $(20, 20)$ and $(40, 40)$), fewer even
than with the exact inverse of $A_\delta$ ($11$ and $10$), but it does
not converge within $100$ iterations at $(80, 80)$, so it is not
robust.""",
r"""A monolithic cycle in which only the patch factorisations and the
coarser levels use $A_\delta$, while the finest-level smoother and
residual use the unshifted operator $A_0$, is a different
preconditioner.  It is what PETSc's multigrid does by default when it
is given $A_0$ as the operator and $A_\delta$ only as the
preconditioning matrix, as it is inside Firedrake's auxiliary-operator
preconditioner; the option \opt{pc_use_amat:false} makes it use
$A_\delta$ throughout, as in Table~\ref{tab:monolithic}.  That cycle
is a multigrid method for $A_0$ rather than an approximation of
$A_\delta^{-1}$.  It needs fewer iterations on the two smaller meshes
($6$ and $8$ at $(20, 20)$ and $(40, 40)$), fewer even than the exact
inverse of $A_\delta$ ($11$ and $10$), but it does not converge within
$100$ iterations at $(80, 80)$, so it is not robust.""")
add("R7 Table 5 caption",
r"""    cycle per application in both (Chebyshev(6) smoothing with the
    column-star patch and RCM-reordered patch factorisations on all
    three levels).  FGMRES with right preconditioning to a relative
    residual of $10^{-8}$; aspect ratio $\varepsilon = 1.3\times10^{-2}$,
    $\Delta t = 100\,\mathrm{s}$, $\delta = 10^{-4}$, $16$ MPI processes.
    $N_x$ is the number of columns of the coarsest multigrid level (the
    finest has $4N_x$); patch DOFs are those of the finest level; times
    are medians of three runs.}""",
r"""    cycle per application in both (three levels: Chebyshev(6) smoothing
    with the column-star patch and RCM-reordered patch factorisations on
    the two finer levels, a direct solve on the coarsest).  FGMRES with
    right preconditioning to a relative residual of $10^{-8}$; aspect
    ratio $\varepsilon = 1.3\times10^{-2}$, $\Delta t = 100\,\mathrm{s}$,
    $\delta = 10^{-4}$ ($C_1 = \delta\,\Delta t = 10^{-2}$), $16$ MPI
    processes.  $N_x$ is the number of columns of the coarsest multigrid
    level (the finest has $4N_x$); the patches have the same size on
    both smoothed levels; times are medians of three runs.}""")
add("R9 conclusion sentence (and stale comment)",
r"""monolithic column-patch smoother applied directly to $A_\delta$
converges in the same number of iterations as the Schur reduction, as
it must, since eliminating the pressure from its patches gives the
fieldsplit patches; the fieldsplit patches are about $20\%$ smaller,
so the Schur reduction is the cheaper of the two.
%  so that the cumulative
% preconditioner cost is lower for the monolithic solver at every
% resolution tested.
""",
r"""monolithic column-patch smoother applied directly to $A_\delta$
returns the same velocity as the Schur reduction in each application,
since eliminating the pressure from its patches gives the fieldsplit
patches, and converges in the same number of iterations to within one;
the fieldsplit patches are about $20\%$ smaller, and the Schur
reduction is the cheaper of the two, by $10$--$15\%$ in preconditioner
time at four of the five resolutions tested.
""")
add("R10 Appendix B relaxation paragraph",
r"""with the damping factor $\theta$.  Neither has been
tuned to be optimised, which is the sense in which the comparison of
Section~\ref{sec:timing} between the fieldsplit and monolithic
smoothers is untuned on both sides.  The coarse grid is solved
directly.""",
r"""with the damping factor $\theta$.  The timings of
Section~\ref{sec:timing} use \opt{ksp_type:chebyshev} with
\opt{ksp_max_it:6} instead, six Chebyshev steps on the interval
estimated by PETSc, for the fieldsplit and the monolithic smoother
alike.  Neither choice has been tuned.  The coarse grid is solved
directly.""")
add("R11 Section 5.2: mesh convention",
r"""LU factorisations for better ordering of the sparse matrix.
""",
r"""LU factorisations for better ordering of the sparse matrix.
The multigrid hierarchy has three levels: a base mesh of $N_x$
intervals, refined twice in the horizontal direction, with $N_z$ layers
on every level.  The mesh sizes $(N_x, N_z)$ quoted in this section
refer to the coarsest level, so that the finest mesh, on which the
system is solved, has $4N_x \times N_z$ cells and
$\Delta x = L_x/(4N_x)$.
""")
add("R12 introduction: 'high cost' and the citation",
r"""This patch structure differs from the monolithic 
patch smoothers applied directly to the full saddle-point system~\cite{farrellMitchellWechsung2019}, in which each patch contains velocity and pressure degrees of freedom. 
Monolithic patches are robust but carry a high cost per application because of their size. """,
r"""This patch structure differs from the monolithic 
patch smoothers applied directly to the full saddle-point system~\cite{PCPATCH}, in which each patch contains velocity and pressure degrees of freedom. 
Monolithic patches are robust but larger, and therefore more expensive to apply. """)
add("R13 Section 5.6: FMW2019 attribution",
r"""monolithic construction is the column-patch instance of the approach
of Farrell, Mitchell, and Wechsung~\cite{farrellMitchellWechsung2019}
discussed in Section~\ref{sec:inner-solvers}.""",
r"""monolithic construction is a column-patch instance of the monolithic
patch relaxation of the PCPATCH framework~\cite{PCPATCH} of
Section~\ref{sec:inner-solvers}.""")
add("R14 Section 3.4: which vertical facets",
r"""those columns together with the horizontal-velocity degrees of freedom
on their vertical facets, problem""",
r"""those columns together with the horizontal-velocity degrees of freedom
on the vertical facet they share, problem""")
add("R15 RCM paragraph: per-layer cost is not constant (optional, pre-existing)",
r"""bandwidth is independent of $N_z$ and the PCApply cost per layer,
reported in the right-most column of each block, is essentially
constant, consistent with a linear-in-$N_z$ scaling of the total cost
(a least-squares fit of PCApply against $N_z$""",
r"""bandwidth is independent of $N_z$ and the PCApply cost grows
linearly in $N_z$: the cost per layer, reported in the right-most
column of each block, falls from $0.169\,\mathrm{s}$ to
$0.100\,\mathrm{s}$ as a fixed overhead is amortised
(a least-squares fit of PCApply against $N_z$""")
add("R16 conclusion: per-layer cost (optional, pre-existing)",
r"""grows superlinearly, while with it the cost per layer is constant.  And a""",
r"""grows superlinearly, while with it the cost grows only linearly.  And a""")
add("R17 Section 5.2: 55pt overfull box (optional, pre-existing)",
r"""smoothing via Firedrake's \texttt{ASMStarPC} 
(\texttt{pc\_star\_construct\_dim:0}).""",
r"""smoothing via Firedrake's \opt{ASMStarPC} 
(\opt{pc_star_construct_dim:0}).""")

add("R18 RCM paragraph: Table 4 runs used GMRES (optional)",
r"""reused across all FGMRES iterations and all time steps.""",
r"""reused across all Krylov iterations and all time steps.""")
add("R19 Section 3.4: smoother actually used in the timings (optional)",
r"""the corrections from all columns are added with a Richardson damping
factor of $\theta = 0.2$. """,
r"""the corrections from all columns are added with a Richardson damping
factor of $\theta = 0.2$ (the timings of Section~\ref{sec:timing} use a
Chebyshev iteration instead; Appendix~\ref{app:petsc-options}). """)

log = []
for name, old, new in R:
    n = t.count(old)
    log.append({"name": name, "count_in_source": n, "old": old, "new": new})
    if n != 1:
        print("MISMATCH", name, n); continue
    t = t.replace(old, new)
# R20 layout: move the Table-5 float before the paragraph that discusses it
b = t.index("\\begin{table}[t]\n    \\centering\n    \\caption{Fieldsplit Schur-complement preconditioning versus a")
e = t.index("\\end{table}\n", b) + len("\\end{table}\n")
block = t[b:e]; t = t[:b] + t[e:]
i = t.index("Table~\\ref{tab:monolithic} reports, for the geometry, time step and")
t = t[:i] + block + "\n" + t[i:]
log.append({"name": "R20 move the Table-5 float block to just before the paragraph 'Table~\\ref{tab:monolithic} reports, ...'", "count_in_source": 1})
open(dst, "w").write(t)
json.dump(log, open(out_json, "w"), indent=1)
print("applied", sum(1 for l in log if l["count_in_source"] == 1), "of", len(log))
