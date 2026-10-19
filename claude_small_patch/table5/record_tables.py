'''
Markdown result tables for claude_paper_update_and_table5.md, generated from the run folders (via parse_table5.py),
so that no number is copied by hand. Prints the tables; usage: python record_tables.py > results_record.md
Counts are "step 1 / warm-up / official" (the three solves of a --tmax 200 --timing run); PCApply is the
Official-Run stage of -log_view (seconds); a "-" means the run is missing or did not finish.
'''
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import parse_table5 as pt          # noqa: E402

R = {r['tag']: r for r in pt.parse_runs()}
OLD = {r['tag']: r for r in pt.parse_old()}
SMALL = [(10, 10), (20, 20), (40, 40), (40, 80)]


def c(tag):
    r = R.get(tag)
    if not r or not r['its_all']:
        return '-'
    s = r['its_all']
    if r['reason'] and r['reason'] != 'CONVERGED_RTOL':
        s += ' (cap)' if r['reason'] == 'CONVERGED_ITS' else f" ({r['reason']})"
    return s


def t(tag, nd=2):
    r = R.get(tag)
    return '-' if not r or r['pcapply_s'] is None else f"{r['pcapply_s']:.{nd}f}"


def table(head, rows):
    out = ['| ' + ' | '.join(head) + ' |', '|' + '---|' * len(head)]
    out += ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows]
    return '\n'.join(out)


def main():
    P = []
    P.append('#### (a) Reproduction of Table 5, serial and parallel (Set A; GMRES left, the 2026-06-01 configuration)\n')
    rows = []
    for nx, nz in SMALL:
        o_fs, o_mo = OLD.get(f'old_fs_nz{nz}'), OLD.get(f'old_mono_nz{nz}')
        rows.append([f'({nx},{nz})', 'FS', o_fs['its'] if o_fs else '-'] +
                    [c(f'A_fs_np{p}_nx{nx}_nz{nz}') for p in (1, 2, 16)] +
                    [t(f'A_fs_np{p}_nx{nx}_nz{nz}') for p in (1, 2, 16)] + [f"{o_fs['pcapply_s']:.2f}" if o_fs else '-'])
        rows.append(['', 'MONO (as run)', o_mo['its'] if o_mo else '-'] +
                    [c(f'A_mono_np{p}_nx{nx}_nz{nz}') for p in (1, 2, 16)] +
                    [t(f'A_mono_np{p}_nx{nx}_nz{nz}') for p in (1, 2, 16)] + [f"{o_mo['pcapply_s']:.2f}" if o_mo else '-'])
    P.append(table(['(N_x,N_z)', 'solver', 'its 2026-06-01 log', 'its serial', 'its 2 ranks', 'its 16 ranks',
                    'PCApply s serial', '2 ranks', '16 ranks', '2026-06-01 (2 ranks)'], rows))
    P.append('\n(10,10) cannot run on 16 ranks: the coarse base mesh has 10 cells (MPI error in the mesh distribution).\n')

    P.append('#### (b) Norm of the stopping test (Set A\', 2 ranks)\n')
    rows = []
    for nx, nz in SMALL:
        rows.append([f'({nx},{nz})', c(f'A_fs_np2_nx{nx}_nz{nz}'), c(f'Ap_fs_right_np2_nx{nx}_nz{nz}'),
                     c(f'A_mono_np2_nx{nx}_nz{nz}'), c(f'Ap_monocheb_right_np2_nx{nx}_nz{nz}'),
                     c(f'Ap_monogmres_right_np2_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'FS, GMRES left', 'FS, GMRES right', 'MONO (as run), GMRES left', 'MONO (as run), FGMRES right',
                    'MONO (as run), GMRES(6) smoother, FGMRES right'], rows))
    P.append('\nThe planned left-preconditioned runs with `-ksp_norm_type unpreconditioned` (`Ap_monocheb_left_rich_*`) '
             'stopped at once: PETSc\'s GMRES does not support the unpreconditioned norm with left preconditioning.\n')

    P.append('#### (c) Converged inner solves (Set B, 2 ranks, outer FGMRES right)\n')
    rows = []
    for nx, nz in [(20, 20), (40, 40)]:
        for tol in ('1e-6', '1e-10'):
            rows.append([f'({nx},{nz})', tol, c(f'B_fs_inner{tol}_np2_nx{nx}_nz{nz}'), c(f'B_mono_inner{tol}_np2_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'inner rtol', 'FS', 'MONO'], rows))
    P.append('\nThe inner FGMRES of MONO is a PCKSP inside `MonolithicShiftPC`; PCKSP uses Pmat by default, so these MONO runs '
             'are converged solves with A_delta (the useAmat leak of Section 2.3 does not occur here).\n')

    P.append('#### (d) Accuracy of the monolithic patch LU (Set Y, 2 ranks, GMRES left)\n')
    rows = []
    for nx, nz in [(20, 20), (40, 40)]:
        rows.append([f'({nx},{nz})', c(f'A_mono_np2_nx{nx}_nz{nz}'), c(f'Y_mono_natural_np2_nx{nx}_nz{nz}'),
                     c(f'Y_mono_rcm_mumps_np2_nx{nx}_nz{nz}'), c(f'Ap_monocheb_right_np2_nx{nx}_nz{nz}'),
                     c(f'Y_mono_natural_right_np2_nx{nx}_nz{nz}'), c(f'A_fs_np2_nx{nx}_nz{nz}'), c(f'Y_fs_rcm_mumps_np2_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'MONO RCM + PETSc LU', 'MONO natural + PETSc LU', 'MONO MUMPS', 'MONO RCM, FGMRES right',
                    'MONO natural, FGMRES right', 'FS RCM + PETSc LU', 'FS MUMPS'], rows))
    P.append('\n(MONO as run in all columns; GMRES left unless stated. Timings of Set Y are not comparable: it ran next to '
             'other jobs; there the MUMPS patch solves were 3-7 x slower than the banded PETSc LU with RCM. '
             '`Y_mono_rcm_mumps_np2_nx80_nz80` was stopped by hand during its first solve (it would have taken about 100 min '
             'on 2 ranks); replaced by the 16-rank run of (g).)\n')

    P.append('#### (e) Larger meshes, 2026-06-01 configuration (Set C, GMRES left)\n')
    rows = []
    for nx, nz in [(80, 80), (80, 160)]:
        rows.append([f'({nx},{nz})', c(f'C_fs_np2_nx{nx}_nz{nz}'), c(f'C_fs_np16_nx{nx}_nz{nz}'),
                     c(f'C_mono_np2_nx{nx}_nz{nz}'), c(f'C_mono_np16_nx{nx}_nz{nz}'),
                     t(f'C_fs_np16_nx{nx}_nz{nz}'), t(f'C_mono_np16_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'FS 2 ranks', 'FS 16 ranks', 'MONO (as run) 2 ranks', 'MONO (as run) 16 ranks',
                    'FS PCApply 16 ranks', 'MONO PCApply 16 ranks'], rows))
    P.append('\nThe rest of Set C (FS and MONO at (80,160) on 2 ranks, MONO at (80,160) on 16 ranks, both solvers at (80,80) '
             'in serial) was dropped once the cause was known: the MONO runs that hit the cap took 8 min (16 ranks) and '
             '20 min (2 ranks) at (80,80), and longer at (80,160) and in serial; the serial question is answered by (a) and '
             '(f). (FS at (80,160) on 16 ranks overlapped two other runs; its PCApply, 5.05 s per application, is 13 % '
             'above the clean T run.)\n')

    P.append('#### (f) The corrected monolithic solver (Set F; `-firedrake_2_monoshiftpc_pc_use_amat false`)\n')
    rows = []
    for nx, nz in SMALL:
        rows.append([f'({nx},{nz})', c(f'A_fs_np2_nx{nx}_nz{nz}'), c(f'F_fspmat_np2_nx{nx}_nz{nz}'),
                     c(f'F_monofix_np2_nx{nx}_nz{nz}'), c(f'Ap_fs_right_np2_nx{nx}_nz{nz}'),
                     c(f'F_monofix_right_np2_nx{nx}_nz{nz}'), c(f'F_monofix_np1_nx{nx}_nz{nz}'),
                     c(f'F_monofix_np16_nx{nx}_nz{nz}'), c(f'A_mono_np2_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'FS left', 'FS (pmat) left', 'MONO fixed left', 'FS right (GMRES)', 'MONO fixed right (FGMRES)',
                    'MONO fixed serial', 'MONO fixed 16 ranks', 'MONO as run left'], rows))
    rows = []
    for nx, nz in [(20, 20), (40, 40)]:
        rows.append([f'({nx},{nz})', c(f'F_fsdirect_np1_nx{nx}_nz{nz}'), c(f'B_fs_inner1e-10_np2_nx{nx}_nz{nz}'),
                     c(f'F_monofixgmres_right_np2_nx{nx}_nz{nz}'), c(f'Ap_monogmres_right_np2_nx{nx}_nz{nz}')])
    P.append('\n' + table(['(N_x,N_z)', 'FS with exact LU on S (serial, GMRES left)', 'FS inner 1e-10 (FGMRES right)',
                           'MONO fixed, GMRES(6) smoother, FGMRES right', 'MONO as run, same'], rows) + '\n')

    P.append('#### (g) The 2026-06-01 monolithic configuration on the new meshes (Set H, 16 ranks)\n')
    rows = []
    for nx, nz in [(20, 20), (40, 40), (80, 80)]:
        rows.append([f'({nx},{nz})', c(f'H_mono_right_np16_nx{nx}_nz{nz}'), c(f'T_monofix_rep1_np16_nx{nx}_nz{nz}'),
                     c(f'T_fs_rep1_np16_nx{nx}_nz{nz}')])
    P.append(table(['(N_x,N_z)', 'MONO as run, FGMRES right', 'MONO fixed, FGMRES right', 'FS, FGMRES right'], rows))
    P.append(f"\nMONO as run at (80,80), 16 ranks: with MUMPS patch LU (GMRES left) {c('H_mono_mumps_np16_nx80_nz80')}, as with "
             f"PETSc's LU (17/100/100 in (e)); with a residual monitor (FGMRES right) {c('H_mono_monitor_right_np16_nx80_nz80')}: "
             "the true residual falls fast (relative 1e-6 after 10-16 iterations) and then stagnates, at 2.0e-8 after 100 "
             "iterations in the first solve and 4.7e-7 in the other two, above the tolerance 1e-8.\n")
    print('\n'.join(P))


if __name__ == '__main__':
    main()
