r'''Proposed old -> new edits for ex_article.tex (sections sec:robustness / sec:nonlinear), built from the file's own
lines so that every old string is exact (trailing spaces included). Verifies that each old string occurs exactly
once, applies all edits to a COPY in this folder, writes edits.json and the unified diff. Never touches the paper.'''
import json, os, difflib
HERE = os.path.dirname(os.path.abspath(__file__))
TEX = '/home/mh2520/molin/phd/yr2/Patch-Preconditioner-for-Boussinesq/ex_article.tex'
src = open(TEX).read()
L = src.split('\n')
def lines(a, b):            # 1-based inclusive, exact text
    return '\n'.join(L[a - 1:b])

E = []   # (id, group, line, old, new)

# ---------------- B: statements the (unchanged) data contradict ----------------
E.append(('B1-dz-text', 'B', 3057, lines(3057, 3059),
'As in the horizontal sweep, every run takes the same number of\n'
'iterations, reflecting the fact that the column-patch\n'
'smoother resolves the vertical coupling exactly within each patch.'))

E.append(('B2-dz-caption', 'B', 3067, lines(3067, 3070),
'    held at their nominal values. As in\n'
'    Figure~\\ref{fig:robust-dx}, every curve ends at the same\n'
'    iteration, reflecting the fact that the column-patch smoother\n'
'    resolves the vertical coupling exactly within each patch.}'))

E.append(('B3-shift-cap-text', 'B', 3226, lines(3226, 3228) + '\ngrows;',
'iteration count diverges (the $C_1 = 1$ curve reaches the iteration cap\n'
'for every $\\Delta t \\le 50\\,\\mathrm{s}$, and the $C_1 = 10^{-1}$ curve for\n'
'every $\\Delta t \\le 10\\,\\mathrm{s}$), while for $C_1 \\le 10^{-3}$, and for\n'
'$C_1 = 10^{-2}$ when $\\Delta t \\ge 5\\,\\mathrm{s}$, it remains bounded and\n'
'decreases as $\\Delta t$\n'
'grows;'))

E.append(('B4-shift-small-dt', 'B', 3258, 'counts are large for every $C_1$ and ordered monotonically in it.',
'counts exceed their values at $\\Delta t \\ge 50\\,\\mathrm{s}$ for every\n'
'$C_1$ and increase monotonically with it.'))

E.append(('B5-shift-caption', 'B', 3298, lines(3298, 3300),
'    near $10^{-2}$, for $\\Delta t \\gg 2/N$.  The curves separate\n'
'    accordingly into those that respect it, which level off beyond\n'
'    $\\Delta t \\sim 2/N$, and those that do not.  Runs that reach the\n'
'    iteration cap of $100$ are omitted.}'))

E.append(('B6-snes-C1', 'B', 3330, 'preconditioned by $A_\\delta^{-1}$ as before. The convergence',
'preconditioned by $A_\\delta^{-1}$ as before, with\n'
'$\\delta = C_1/\\Delta t$ and $C_1 = 10^{-5}$ in the resolution and\n'
'aspect-ratio sweeps of Figures~\\ref{fig:snes-dx}--\\ref{fig:snes-ar}. The convergence'))

E.append(('B7-snes-ar-text', 'B', 3367, lines(3367, 3369),
'$\\Delta t = 5000\\,\\mathrm{s}$. At $\\Delta t = 5000\\,\\mathrm{s}$ the\n'
'curves are tightly clustered; at $\\Delta t = 100\\,\\mathrm{s}$ the\n'
'cumulative count falls from $19$--$20$ at $L_z = 12$ and $8\\,\\mathrm{km}$\n'
'to $11$ at $L_z = 1\\,\\mathrm{km}$, the same improvement in thinner\n'
'domains as in the linear sweep (Remark~\\ref{rem:thin-faster}).  In\n'
'neither case does the count grow as $\\varepsilon$ decreases, confirming\n'
'that the preconditioner retains its aspect-ratio robustness in the\n'
'nonlinear setting.'))

E.append(('B8-snes-ar-caption', 'B', 3384, lines(3384, 3386),
'    (left) and $\\Delta t = 5000\\,\\mathrm{s}$ (right) for five\n'
'    aspect ratios, $L_z \\in \\{12, 8, 4, 2, 1\\}\\,\\mathrm{km}$\n'
'    ($\\varepsilon = L_z/L_x$ from $4\\times10^{-2}$ to $3.3\\times10^{-3}$).\n'
'    Thinner domains converge in fewer iterations at\n'
'    $\\Delta t = 100\\,\\mathrm{s}$ ($11$ against $19$--$20$ cumulative\n'
'    iterations); at $\\Delta t = 5000\\,\\mathrm{s}$ the curves are tightly\n'
'    clustered.  In the right panel the final Newton step, which takes\n'
'    the error below $10^{-12}$, is not plotted.}'))

E.append(('B9-snes-shift-text', 'B', 3397, lines(3397, 3398),
'of Section~\\ref{sec:robustness}: for $\\Delta t \\le 500\\,\\mathrm{s}$ the\n'
'iteration count decreases monotonically as $C_1$ decreases, so that the\n'
'smallest value tested, $C_1 = 10^{-6}$, gives the lowest counts, while\n'
'for $C_1 \\ge 5\\times10^{-4}$ and $\\Delta t \\le 100\\,\\mathrm{s}$ the count\n'
'exceeds $100$; for $\\Delta t \\ge 10^{4}\\,\\mathrm{s}$ all seven values give\n'
'nearly the same count; and'))

E.append(('B10-snes-shift-hist-caption', 'B', 3410, lines(3410, 3412),
'        \\caption{Nonlinear shift-parameter sensitivity: convergence\n'
'        history at $\\Delta t = 500\\,\\mathrm{s}$ for the seven values\n'
'        $C_1 \\in \\{10^{-3}, 5\\times10^{-4}, 10^{-4}, 5\\times10^{-5},\n'
'        10^{-5}, 5\\times10^{-6}, 10^{-6}\\}$.}'))

E.append(('B11-snes-shift-sweep-caption', 'B', 3418, lines(3418, 3420),
'        \\caption{Nonlinear shift-parameter sweep: total KSP\n'
'        iteration count over the Newton steps vs $\\Delta t$ for the\n'
'        seven values of $C_1$ of Figure~\\ref{fig:snes-shift-hist}.\n'
'        Runs in which a Newton step reached the cap of $100$ KSP\n'
'        iterations, or which did not complete, are omitted.}'))

# ---------------- C: count-convention harmonisation (go with B11) ----------------
E.append(('C1-fig12-true-totals', 'C', 3417,
'        \\includegraphics[width=\\linewidth]{Image_restyled/error_shift_nonlinear}',
'        \\includegraphics[width=\\linewidth]{Image_restyled/corrected/error_shift_nonlinear}'))

E.append(('C2-linear-history-note', 'C', 3050, lines(3050, 3051),
'    fixed $\\Delta t = 50\\,\\mathrm{s}$, with remaining parameters\n'
'    held at their nominal values.  Here and in\n'
'    Figures~\\ref{fig:robust-dz} and~\\ref{fig:robust-ar} the converged\n'
'    iterate (error below $10^{-13}$) is not plotted, so the curve of a\n'
'    run with $n$ iterations ends at index $n-1$.}'))

# ---------------- A: only if the legend-fixed candidate figures are adopted ----------------
E.append(('A1-ar-caption-if-short-legend', 'A', 3189, 'for seven aspect ratios at fixed\n    $\\Delta t = 50\\,\\mathrm{s}$.',
'for the seven depths\n'
'    $L_z \\in \\{12, 8, 4, 2, 1, 0.5, 0.25\\}\\,\\mathrm{km}$ of\n'
'    Table~\\ref{tab:mode-damping} at fixed $\\Delta t = 50\\,\\mathrm{s}$.'))

E.append(('A2-fig11-caption-if-no-legend', 'A', 3412, None, None))   # handled as a variant of B10, see report

out = []
copy = src
for eid, grp, ln, old, new in E:
    if old is None:
        continue
    n = src.count(old); m = copy.count(old)
    status = 'OK' if (n == 1 and m == 1) else f'COUNT src={n} copy={m}'
    print(f'{eid:32s} line {ln:5d} occurrences {n} {status}')
    assert n == 1 and m == 1, eid
    copy = copy.replace(old, new, 1)
    out.append(dict(id=eid, group=grp, line=ln, old=old, new=new))
json.dump(out, open(os.path.join(HERE, 'edits.json'), 'w'), indent=1)
open(os.path.join(HERE, 'ex_article.figure_referee.tex'), 'w').write(copy)
diff = difflib.unified_diff(src.split('\n'), copy.split('\n'), 'ex_article.tex', 'ex_article.figure_referee.tex', n=0, lineterm='')
open(os.path.join(HERE, 'edits.diff'), 'w').write('\n'.join(diff) + '\n')
# brace balance check of the changed text
for d in out:
    for k in ('old', 'new'):
        s = d[k].replace('\\{', '').replace('\\}', '')
        if s.count('{') - s.count('}') != (d['old'].replace('\\{', '').replace('\\}', '').count('{') - d['old'].replace('\\{', '').replace('\\}', '').count('}')):
            print('BRACE BALANCE CHANGED in', d['id'], k)
print('total braces old/new file:', src.count('{') - src.count('}'), copy.count('{') - copy.count('}'))
