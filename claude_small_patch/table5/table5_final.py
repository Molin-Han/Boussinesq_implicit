'''
The new Table 5 from the timed runs of run_table5_fixed.sh (Set T: 16 ranks, FGMRES right, three repetitions).

Per solver (fs, fspmat, monofix) and mesh: outer iterations of the official solve (checked equal across the
repetitions), Official-Run PCApply time (median and range of the repetitions), time per iteration, finest-level patch
size, total number of unknowns. Writes results_final.md and table5_tabular.tex (the tabular for the paper).
Usage: python table5_final.py
'''
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import parse_table5 as pt          # noqa: E402

MESHES = [(20, 20), (40, 40), (80, 80), (80, 160), (160, 160)]
SOLVERS = [('fs', 'fieldsplit Schur (Amat = MatSchurComplement, the paper\'s configuration)'),
           ('fspmat', 'fieldsplit Schur, assembled auxiliary operator on every level (pc_use_amat false)'),
           ('monofix', 'monolithic on A_delta on every level (pc_use_amat false)')]


def ndofs(nx, nz):
    n = 4 * nx                     # two horizontal refinements of the coarse base mesh
    return 20 * n * nz + 4 * n     # RT_2 (u_xz) + DG1 (u_y) + W_theta,2 (b) + DG1 (p), periodic, degree 2


def collect():
    rows = {r['tag']: r for r in pt.parse_runs()}
    out = {}
    for s, _ in SOLVERS:
        for (nx, nz) in MESHES:
            reps = [rows.get(f'T_{s}_rep{k}_np16_nx{nx}_nz{nz}') for k in (1, 2, 3)]
            reps = [r for r in reps if r and r['its'] is not None and r['pcapply_s'] is not None]
            if not reps:
                continue
            its = sorted({r['its'] for r in reps})
            t = [r['pcapply_s'] for r in reps]
            k = [r['ksp_s'] for r in reps]
            out[(s, nx, nz)] = dict(its=its, n=len(reps), t_med=st.median(t), t_min=min(t), t_max=max(t),
                                    ksp_med=st.median(k), patch=max(r['patch'] for r in reps),
                                    reasons=sorted({r['reason'] for r in reps}), its_all=[r['its_all'] for r in reps])
    hyb = {}
    for (nx, nz) in [(20, 20), (40, 40), (80, 80)]:
        r = rows.get(f'H_mono_right_np16_nx{nx}_nz{nz}')
        if r:
            hyb[(nx, nz)] = r
    return out, hyb


def main():
    out, hyb = collect()
    L = ['# New Table 5: 16 MPI ranks, FGMRES (right preconditioning), rtol 1e-8, median of the repetitions', '',
         '| (N_x, N_z) | fine mesh | unknowns | solver | patch dofs | its (official) | PCApply s median [min, max] (n) | s per it | KSPSolve s |',
         '|---|---|---|---|---|---|---|---|---|']
    for (nx, nz) in MESHES:
        for s, _ in SOLVERS:
            d = out.get((s, nx, nz))
            if not d:
                continue
            its = d['its'][0] if len(d['its']) == 1 else '/'.join(map(str, d['its']))
            per = d['t_med'] / d['its'][0]
            L.append(f"| ({nx},{nz}) | {4 * nx} x {nz} | {ndofs(nx, nz)} | {s} | {d['patch']} | {its} | "
                     f"{d['t_med']:.2f} [{d['t_min']:.2f}, {d['t_max']:.2f}] ({d['n']}) | {per:.3f} | {d['ksp_med']:.2f} |")
    L += ['', 'All three solves of each run (step 1 / warm-up / official):', '']
    for key, d in sorted(out.items()):
        L.append(f"* {key}: {', '.join(d['its_all'])}; reasons {d['reasons']}")
    if hyb:
        L += ['', '## The 2026-06-01 monolithic configuration (finest-level smoother on the unshifted Jacobian), FGMRES right, 16 ranks', '',
              '| (N_x, N_z) | its (3 solves) | PCApply s |', '|---|---|---|']
        for (nx, nz), r in sorted(hyb.items()):
            L.append(f"| ({nx},{nz}) | {r['its_all']} | {r['pcapply_s'] if r['pcapply_s'] is None else round(r['pcapply_s'], 2)} |")
    open(os.path.join(HERE, 'results_final.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
    # LaTeX tabular: paper configuration (fs) against the corrected monolithic solver (monofix)
    T = [r'    \small', r'    \begin{tabular}{cc ccc ccc}', r'        \toprule',
         r'        & & \multicolumn{3}{c}{Fieldsplit Schur} & \multicolumn{3}{c}{Monolithic} \\',
         r'        \cmidrule(lr){3-5} \cmidrule(lr){6-8}',
         r'        $(N_x, N_z)$ & Unknowns & Patch DOFs & PCApply & Its',
         r'                     & Patch DOFs & PCApply & Its \\', r'        \midrule']
    for (nx, nz) in MESHES:
        a, b = out.get(('fs', nx, nz)), out.get(('monofix', nx, nz))
        if not (a and b):
            continue
        n = ndofs(nx, nz)
        ns = f'{n / 1e6:.2f}M' if n >= 1e6 else f'{n / 1e3:.0f}k'
        T.append(f"        $({nx}, {nz})$ & {ns} & {a['patch']} & {a['t_med']:.2f}s & {a['its'][0]} & "
                 f"{b['patch']} & {b['t_med']:.2f}s & {b['its'][0]} \\\\")
    T += [r'        \bottomrule', r'    \end{tabular}']
    open(os.path.join(HERE, 'table5_tabular.tex'), 'w').write('\n'.join(T) + '\n')
    print('\n'.join(T))


if __name__ == '__main__':
    main()
