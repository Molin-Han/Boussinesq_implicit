'''
Parse the Table-5 re-runs in runs/<tag>/ (see run_table5.sh) and the eight 2026-06-01 logs behind the paper's
table, and write results.csv + results.md.

Per run: outer iteration counts of the three solves of the run (step 1, warm-up, official; last = official),
Official-Run stage PCApply count and time and KSPSolve time from -log_view, the largest patch size printed by
pc_*_view_patch_sizes, the wall time of the whole run and the process count from the log header.
Usage: python parse_table5.py [--old]   (--old: also tabulate yr2/data_firedrake26_{timing,mono}/log_view_*.txt)
'''
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'runs')
OLD = os.path.join(HERE, '..', '..', '..')          # yr2/
TAG = re.compile(r'^(?P<set>[A-Z][a-z]?)_(?P<solver>[a-z]+)(?P<variant>(?:_[A-Za-z0-9.+-]+?)*?)'
                 r'_np(?P<np>\d+)_nx(?P<nx>\d+)_nz(?P<nz>\d+)$')


def stage_events(logfile, stage='Official-Run', events=('PCApply', 'KSPSolve', 'MatMult', 'PCApplyOnBlocks')):
    out = {}
    if not os.path.exists(logfile):
        return out
    inside = False
    with open(logfile) as f:
        for line in f:
            if line.startswith('--- Event Stage'):
                inside = stage in line
                continue
            if inside:
                parts = line.split()
                if len(parts) > 4 and parts[0] in events and parts[0] not in out:
                    out[parts[0]] = (int(parts[1]), float(parts[3]))
    return out


def header_nprocs(logfile):
    if not os.path.exists(logfile):
        return None
    with open(logfile) as f:
        for line in f:
            m = re.search(r'with (\d+) process', line)
            if m:
                return int(m.group(1))
    return None


def parse_runlog(runlog):
    '''outer counts of the three solves, largest patch size, reasons, and the inner iteration counts
    (fieldsplit_1 or monoshiftpc_ksp converged lines) grouped by outer solve'''
    its, patch, reasons, inner, cur = [], 0, [], [], []
    if not os.path.exists(runlog):
        return its, patch, reasons, inner
    with open(runlog) as f:
        for line in f:
            m = re.search(r'Linear (\S+) solve (?:converged|did not converge) due to (\S+) iterations (\d+)', line)
            if m:
                if m.group(1) == 'firedrake_2_':
                    reasons.append(m.group(2)); its.append(int(m.group(3))); inner.append(cur); cur = []
                else:
                    cur.append(int(m.group(3)))
            m = re.search(r'maximum patch sizes : \d+ / [\d.]+ / (\d+)', line)
            if m:
                patch = max(patch, int(m.group(1)))
    return its, patch, reasons, inner


def parse_runs():
    rows = []
    for d in sorted(glob.glob(os.path.join(RUNS, '*'))):
        tag = os.path.basename(d)
        m = TAG.match(tag)
        if not m:
            print('skip', tag, file=sys.stderr); continue
        g = m.groupdict()
        g['variant'] = g['variant'].lstrip('_')
        mi = re.search(r'inner([0-9]+e[-+]?[0-9]+)', g['variant'])
        g['inner'] = mi.group(1) if mi else ''
        its, patch, reasons, inner = parse_runlog(os.path.join(d, 'run.log'))
        if not patch:   # with --timing the patch sizes are printed into the ksp_view file
            patch = parse_runlog(os.path.join(d, 'slice3D.txt'))[1]
        inner_last = inner[-1] if inner else []
        ev = stage_events(os.path.join(d, 'log.txt'))
        st = os.path.join(d, 'status')
        rc, wall = (open(st).read().split() + ['', ''])[:2] if os.path.exists(st) else ('', '')
        rows.append(dict(tag=tag, set=g['set'], solver=g['solver'], variant=g['variant'] or '', inner=g['inner'] or '',
                         np=int(g['np']), nx=int(g['nx']), nz=int(g['nz']), rc=rc, wall_s=wall,
                         its=its[-1] if its else None, its_all='/'.join(map(str, its)), reason=reasons[-1] if reasons else '',
                         patch=patch, pcapply_n=ev.get('PCApply', (None, None))[0], pcapply_s=ev.get('PCApply', (None, None))[1],
                         ksp_s=ev.get('KSPSolve', (None, None))[1], matmult_n=ev.get('MatMult', (None, None))[0],
                         blocks_n=ev.get('PCApplyOnBlocks', (None, None))[0], nprocs_log=header_nprocs(os.path.join(d, 'log.txt')),
                         inner_sum=sum(inner_last) if inner_last else '', inner_n=len(inner_last) if inner_last else '',
                         inner_list='/'.join(map(str, inner_last))))
    return rows


def parse_old():
    rows = []
    for solver, folder in (('fs', 'data_firedrake26_timing'), ('mono', 'data_firedrake26_mono')):
        for nz in (10, 20, 40, 80):
            f = os.path.join(OLD, folder, f'log_view_{nz}.txt')
            ev = stage_events(f)
            n = ev.get('PCApply', (None, None))[0]
            rows.append(dict(tag=f'old_{solver}_nz{nz}', set='old', solver=solver, variant='', inner='', np=header_nprocs(f),
                             nx={10: 10, 20: 20, 40: 40, 80: 40}[nz], nz=nz, rc='', wall_s='', its=(n - 1) if n else None,
                             its_all='(PCApply-1)', reason='', patch='', pcapply_n=n, pcapply_s=ev.get('PCApply', (None, None))[1],
                             ksp_s=ev.get('KSPSolve', (None, None))[1], matmult_n=ev.get('MatMult', (None, None))[0],
                             blocks_n=ev.get('PCApplyOnBlocks', (None, None))[0], nprocs_log=header_nprocs(f),
                             inner_sum='', inner_n='', inner_list=''))
    return rows


def fmt(x, nd=2):
    if x is None or x == '':
        return '-'
    return f'{x:.{nd}f}' if isinstance(x, float) else str(x)


def main():
    rows = parse_runs()
    if '--old' in sys.argv:
        rows = parse_old() + rows
    keys = ['tag', 'set', 'solver', 'variant', 'inner', 'np', 'nx', 'nz', 'its', 'its_all', 'reason', 'patch',
            'pcapply_n', 'pcapply_s', 'ksp_s', 'matmult_n', 'blocks_n', 'wall_s', 'rc', 'nprocs_log', 'inner_n', 'inner_sum', 'inner_list']
    with open(os.path.join(HERE, 'results.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, '') for k in keys})
    lines = ['| tag | np | (Nx,Nz) | its (3 solves) | inner its (official) | patch dofs | PCApply n | PCApply s | s/it | KSPSolve s | wall s |',
             '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        per = (r['pcapply_s'] / r['its']) if (r['pcapply_s'] is not None and r['its']) else None
        lines.append(f"| {r['tag']} | {r['np']} | ({r['nx']},{r['nz']}) | {r['its_all'] or '-'} | {r.get('inner_list') or '-'} | {fmt(r['patch'])} | "
                     f"{fmt(r['pcapply_n'])} | {fmt(r['pcapply_s'])} | {fmt(per, 3)} | {fmt(r['ksp_s'])} | {r['wall_s'] or '-'} |")
    open(os.path.join(HERE, 'results.md'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
