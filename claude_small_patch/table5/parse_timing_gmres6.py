#!/usr/bin/env python3
'''
Summarise the timing runs of run_timing_tables_gmres6.sh (GMRES(6) level smoother, 2026-11-16) for the paper.

Reads <runs>/<tag>/{status,run.log,log.txt} and writes, next to the runs folder:
  results_gmres6.csv   one row per run
  results_gmres6.md    medians over the repetitions, the numbers the text quotes, and draft LaTeX rows for
                       tab:rcm and tab:monolithic in the layout of the paper's current tables

Sets (tags as the run script makes them):
  R_{rcm,norcm}_rep<k>_np2_nx20_nz<Nz>       tab:rcm            linear, PCApply of stage "Official-Run"
  T_{fs,mono}_rep<k>_np16_nx<Nx>_nz<Nz>       tab:monolithic     linear, PCApply of stage "Official-Run"
  B_{fs,mono}_inner1e-10_np2_nx<Nx>_nz<Nz>    "exact solve" sentence below tab:monolithic: outer iterations
  C_{star,vanka}_rep<k>_np2_nx20_nz20         single column's cost per outer iteration (nonlinear solver):
                                              PCApply time / PCApply count of stage "Second_Run"
PCApply includes the numerical refactorisation of the patch and coarse matrices (as in the paper). With FGMRES and
right preconditioning there is one preconditioner application per outer iteration.

    python3 parse_timing_gmres6.py [--runs DIR]        (default: runs_gmres6 next to this script)
Pure Python 3 (no numpy needed).
'''
import argparse
import csv
import glob
import os
import re
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))
TAG = re.compile(r'^(?P<set>[RTBC])_(?P<variant>[a-z]+)(?:_inner(?P<inner>[0-9.e+-]+))?(?:_rep(?P<rep>\d+))?'
                 r'_np(?P<np>\d+)_nx(?P<nx>\d+)_nz(?P<nz>\d+)$')
# Columns of the paper's tables that do not depend on the smoother (from the current tables and the formulas
# unknowns = 20 n N_z + 4 n with n = 4 N_x, fieldsplit patch 30 N_z + 8, monolithic patch 38 N_z + 8):
BANDWIDTH = {'rcm': {10: 48, 20: 48, 40: 48, 80: 48}, 'norcm': {10: 531, 20: 1051, 40: 1571, 80: 2091}}
MESHES_T = [(20, 20), (40, 40), (80, 80), (80, 160), (160, 160)]
NZ_R = [10, 20, 40, 80]


def unknowns(nx, nz):
    n = 4 * nx
    return 20 * n * nz + 4 * n


def stage_events(logfile, stage, events=('PCApply', 'KSPSolve', 'SNESSolve', 'PCSetUp')):
    '''{event: (count, max time over ranks)} in the first summary block of the named -log_view stage'''
    out, inside = {}, False
    if not os.path.exists(logfile):
        return out
    with open(logfile, errors='replace') as f:
        for line in f:
            if line.startswith('--- Event Stage'):
                inside = line.rstrip().endswith(': ' + stage)
                continue
            if inside:
                p = line.split()
                if len(p) > 4 and p[0] in events and p[0] not in out:
                    try:
                        out[p[0]] = (int(p[1]), float(p[3]))
                    except ValueError:
                        pass
    return out


def nprocs(logfile):
    if os.path.exists(logfile):
        with open(logfile, errors='replace') as f:
            for line in f:
                m = re.search(r'with (\d+) MPI process|with (\d+) process', line)
                if m:
                    return int(m.group(1) or m.group(2))
    return None


def outer_counts(runlog, prefix):
    its, reasons = [], []
    if os.path.exists(runlog):
        with open(runlog, errors='replace') as f:
            for line in f:
                m = re.search(rf'Linear {prefix} solve (converged|did not converge) due to (\S+) iterations (\d+)', line)
                if m:
                    reasons.append(m.group(2)); its.append(int(m.group(3)))
    return its, reasons


def expected_tags():
    '''the 64 runs of run_timing_tables_gmres6.sh R T B C'''
    t = [f'R_{v}_rep{k}_np2_nx20_nz{nz}' for k in (1, 2, 3) for nz in NZ_R for v in ('rcm', 'norcm')]
    t += [f'T_{v}_rep{k}_np16_nx{nx}_nz{nz}' for k in (1, 2, 3) for nx, nz in MESHES_T for v in ('fs', 'mono')]
    t += [f'B_{v}_inner1e-10_np2_nx{nx}_nz{nz}' for nx, nz in ((20, 20), (40, 40)) for v in ('fs', 'mono')]
    t += [f'C_{p}_rep{k}_np2_nx20_nz20' for k in (1, 2, 3) for p in ('star', 'vanka')]
    return t


def parse(runs):
    rows = []
    for d in sorted(glob.glob(os.path.join(runs, '*'))):
        if not os.path.isdir(d):
            continue
        tag = os.path.basename(d)
        m = TAG.match(tag)
        if not m:
            continue
        g = m.groupdict()
        nl = g['set'] == 'C'
        st_ = open(os.path.join(d, 'status')).read().split() if os.path.exists(os.path.join(d, 'status')) else []
        rc, wall = (st_ + ['', ''])[:2]
        its, reasons = outer_counts(os.path.join(d, 'run.log'), 'firedrake_1_' if nl else 'firedrake_2_')
        ev = stage_events(os.path.join(d, 'log.txt'), 'Second_Run' if nl else 'Official-Run')
        n_app, t_app = ev.get('PCApply', (None, None))
        rows.append(dict(tag=tag, set=g['set'], variant=g['variant'], inner=g['inner'] or '', rep=int(g['rep'] or 0),
                         np=int(g['np']), nx=int(g['nx']), nz=int(g['nz']), rc=rc, wall_s=wall,
                         its=(n_app if nl else (its[-1] if its else None)),   # nonlinear: total over the Newton steps
                         its_all='/'.join(map(str, its)), reason=(reasons[-1] if reasons else ''),
                         pcapply_n=n_app, pcapply_s=t_app,
                         solve_s=ev.get('SNESSolve' if nl else 'KSPSolve', (None, None))[1],
                         np_log=nprocs(os.path.join(d, 'log.txt'))))
    return rows


def med(xs):
    xs = [x for x in xs if x is not None]
    return st.median(xs) if xs else None


def spread(xs):
    xs = [x for x in xs if x is not None]
    if len(xs) < 2:
        return None
    m = st.median(xs)
    return (max(xs) - min(xs)) / m if m else None


def sig(x, n=4):
    '''n significant digits, as in tab:rcm (1.685s, 35.01s, 139.7s)'''
    if x is None:
        return '--'
    s = f'{x:.{n}g}'
    if '.' not in s and 'e' not in s and len(s.lstrip('-')) < n:
        s += '.' + '0' * (n - len(s.lstrip('-')))
    return s


def group(rows, **kw):
    '''the runs with the given fields that finished successfully and have a PCApply time (the others are listed
    under Problems and left out of every median)'''
    return [r for r in rows if all(r[k] == v for k, v in kw.items()) and r['rc'] == '0' and r['pcapply_s'] is not None]


def vals(g, key):
    return sorted({r[key] for r in g if r[key] is not None})


def affine_r2(xs, ys):
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    b = sxy / sxx
    a = my - b * mx
    ss_res = sum((y - a - b * x) ** 2 for x, y in zip(xs, ys))
    ss_tot = sum((y - my) ** 2 for y in ys)
    return a, b, 1 - ss_res / ss_tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--runs', default=os.path.join(HERE, 'runs_gmres6'))
    a = ap.parse_args()
    rows = parse(a.runs)
    outdir = os.path.dirname(os.path.abspath(a.runs))
    keys = ['tag', 'set', 'variant', 'inner', 'rep', 'np', 'np_log', 'nx', 'nz', 'rc', 'wall_s', 'its', 'its_all',
            'reason', 'pcapply_n', 'pcapply_s', 'solve_s']
    with open(os.path.join(outdir, 'results_gmres6.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})

    L = [f'# Timing runs with the GMRES(6) level smoother: summary', '',
         f'Runs folder: `{os.path.abspath(a.runs)}` ({len(rows)} runs). PCApply = time of the preconditioner '
         'applications of the timed solve (stage "Official-Run"; "Second_Run" for set C), maximum over ranks, '
         'including the numerical refactorisation. Medians over the repetitions; spread = (max - min)/median.', '']

    bad = [r for r in rows if r['rc'] != '0' or r['pcapply_s'] is None
           or (r['set'] != 'C' and (not r['reason'].startswith('CONVERGED') or r['reason'] == 'CONVERGED_ITS'))
           or (r['np_log'] is not None and r['np_log'] != r['np'])]
    missing = sorted(set(expected_tags()) - {r['tag'] for r in rows})
    L += ['## Completeness', '',
          f"{len(rows)} run folders; {len(expected_tags()) - len(missing)} of the {len(expected_tags())} runs of the full "
          f"script are present." + ('' if not missing else f" Missing ({len(missing)}): " + ', '.join(f'`{t}`' for t in missing)), '']
    L += ['## Problems', '']
    L += [f"- `{r['tag']}`: rc={r['rc'] or 'missing (unfinished?)'}, reason={r['reason'] or '-'}, "
          f"PCApply={'missing' if r['pcapply_s'] is None else 'ok'}, ranks in log={r['np_log']}"
          + (' (reached the 100-iteration cap)' if r['reason'] == 'CONVERGED_ITS' else '') for r in bad] or ['- none']
    if bad:
        L += ['', 'Runs that failed or did not finish are left out of the medians below; a rank mismatch means the run '
                  "was not made as intended. Rerun such runs (delete the folder, then the run script with ONLY='^<tag>$') "
                  'before using the numbers.']
    L += ['']

    # ---- set R: tab:rcm -------------------------------------------------------------------------------------
    R = group(rows, set='R')
    if R:
        L += ['## Set R: tab:rcm (fieldsplit solver, 2 ranks)', '',
              '| (N_x, N_z) | RCM: its | PCApply n | PCApply s median [min, max] | /N_z | no RCM: its | PCApply s median [min, max] | /N_z |',
              '|---|---|---|---|---|---|---|---|']
        tex = []
        fit_x, fit_y = [], []
        for nz in NZ_R:
            cells = []
            for v in ('rcm', 'norcm'):
                g = group(R, variant=v, nz=nz)
                t = [r['pcapply_s'] for r in g]
                its = vals(g, 'its')
                napp = vals(g, 'pcapply_n')
                mt = med(t)
                cells.append((its, napp, mt, (min(t) if t else None, max(t) if t else None), (mt / nz) if mt else None))
            (ri, rn, rt, rr, rz), (ni, nn, nt, nr, nzz) = cells
            L.append(f"| (20, {nz}) | {ri} | {rn} | {sig(rt)} [{sig(rr[0])}, {sig(rr[1])}] | {sig(rz, 3)} | {ni} | "
                     f"{sig(nt)} [{sig(nr[0])}, {sig(nr[1])}] | {sig(nzz, 3)} |")
            tex.append(f"        $(20, {nz})$ & {BANDWIDTH['rcm'][nz]} & {sig(rt)}s & {sig(rz, 3)}s & "
                       f"{BANDWIDTH['norcm'][nz]} & {sig(nt)}s & {sig(nzz, 3)}s \\\\")
            if rt is not None:
                fit_x.append(nz); fit_y.append(rt)
        if len(fit_x) >= 3:
            a0, b0, r2 = affine_r2(fit_x, fit_y)
            L += ['', f'Affine fit of the RCM PCApply against N_z: {a0:.3f} s + {b0:.4f} s x N_z, R^2 = {r2:.4f}.']
        L += ['', 'Draft rows for tab:rcm (bandwidths unchanged: they depend only on the patch matrices):', '', '```latex'] + tex + ['```', '']

    # ---- set T: tab:monolithic -----------------------------------------------------------------------------
    T = group(rows, set='T')
    if T:
        L += ['## Set T: tab:monolithic (16 ranks, FGMRES right, GMRES(6) + star patch + RCM)', '',
              '| (N_x, N_z) | unknowns | fs its | fs PCApply s median [min, max] (spread) | fs s/unknown | mono its | mono PCApply s median [min, max] (spread) | mono/fs | mono higher by | fs cheaper by |',
              '|---|---|---|---|---|---|---|---|---|---|']
        tex = []
        for nx, nz in MESHES_T:
            out = {}
            for v in ('fs', 'mono'):
                g = group(T, variant=v, nx=nx, nz=nz)
                t = [r['pcapply_s'] for r in g]
                out[v] = (vals(g, 'its'), med(t),
                          (min(t) if t else None, max(t) if t else None), spread(t), len(g))
            (fi, ft, fr, fsp, fn), (mi, mt, mr, msp, mn) = out['fs'], out['mono']
            nu = unknowns(nx, nz)
            ratio = (mt / ft) if (mt and ft) else None
            L.append(f"| ({nx}, {nz}) | {nu} | {fi} | {sig(ft)} [{sig(fr[0])}, {sig(fr[1])}] ({fsp:.1%}) | "
                     f"{(ft / nu):.2e} | {mi} | {sig(mt)} [{sig(mr[0])}, {sig(mr[1])}] ({msp:.1%}) | {ratio:.3f} | "
                     f"{ratio - 1:.1%} | {1 - 1 / ratio:.1%} |"
                     if (ft and mt and fsp is not None and msp is not None) else
                     f"| ({nx}, {nz}) | {nu} | {fi} | {sig(ft)} ({fn} runs) | - | {mi} | {sig(mt)} ({mn} runs) | - | - | - |")
            k = f'{nu / 1e3:.0f}k' if nu < 1e6 else f'{nu / 1e6:.2f}M'
            fits = '/'.join(map(str, fi)) or '--'
            mits = '/'.join(map(str, mi)) or '--'
            tex.append(f"        $({nx}, {nz})$ & {k} & {30 * nz + 8} & {ft:.2f}s & {fits} & {38 * nz + 8} & {mt:.2f}s & {mits} \\\\"
                       if (ft and mt) else f"        % ({nx}, {nz}): runs missing")
        L += ['', '"mono higher by" (mono/fs - 1) is the comparison of the Cost subsection ("11--17% higher"); "fs cheaper by" '
                  '(1 - fs/mono) is the one of the Conclusions ("about 10--15% cheaper"). s/unknown divides by the exact '
                  'number of unknowns (the October text divided by the rounded 3.2e4).', '',
              'Draft rows for tab:monolithic (unknowns and patch DOFs from the formulas above; check "1.03M"/"2.05M"):',
              '', '```latex'] + tex + ['```', '']

    # ---- set B: exact solve ---------------------------------------------------------------------------------
    B = group(rows, set='B')
    if B:
        L += ['## Set B: the cycle iterated to 1e-10 (2 ranks)', '', '| solver | (N_x, N_z) | outer its (step 1 / warm-up / official) |',
              '|---|---|---|']
        for r in sorted(B, key=lambda r: (r['variant'], r['nx'])):
            L.append(f"| {r['variant']} | ({r['nx']}, {r['nz']}) | {r['its_all']} |")
        L += ['', 'Compare the official counts with the one-cycle counts of set T at (20,20) and (40,40).', '']

    # ---- set C: single column ---------------------------------------------------------------------------------
    C = group(rows, set='C')
    if C:
        L += ['## Set C: cost of one outer iteration, star patch against single column (nonlinear solver, 2 ranks)', '',
              '| patch | runs | applications (= FGMRES its of the timed step) | PCApply s median | s per application median [min, max] |',
              '|---|---|---|---|---|']
        per = {}
        for v in ('star', 'vanka'):
            g = group(C, variant=v)
            pa = [r['pcapply_s'] / r['pcapply_n'] for r in g if r['pcapply_s'] and r['pcapply_n']]
            per[v] = med(pa)
            L.append(f"| {v} | {len(g)} | {vals(g, 'pcapply_n')} | {sig(med([r['pcapply_s'] for r in g]))} | "
                     f"{sig(per[v])} [{sig(min(pa) if pa else None)}, {sig(max(pa) if pa else None)}] |")
        if per.get('star') and per.get('vanka'):
            L += ['', f"One application with the single column costs {per['vanka'] / per['star']:.3f} times one with the star "
                      f"patch, i.e. it is {1 - per['vanka'] / per['star']:.0%} cheaper (paper and supplement say 'about 30%'). "
                      f"The single column is then faster only where the two counts are within a factor of "
                      f"{per['star'] / per['vanka']:.2f}.", '']

    L += ['## All runs', '', '| tag | rc | wall s | ranks (log) | its | reason | PCApply n | PCApply s | solve s |', '|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        L.append(f"| {r['tag']} | {r['rc']} | {r['wall_s']} | {r['np_log']} | {r['its_all'] if r['set'] != 'C' else r['its']} | "
                 f"{r['reason']} | {r['pcapply_n']} | {sig(r['pcapply_s'])} | {sig(r['solve_s'])} |")
    open(os.path.join(outdir, 'results_gmres6.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))


if __name__ == '__main__':
    main()
