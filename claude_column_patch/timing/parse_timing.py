'''
Collect the timing runs of run_timing.sh into timing_summary.csv / timing_summary.md.

For every run directory <name>_np<ranks>_<patch>_rep<r>/ it reads
  log_view.txt : PETSc -log_view; the time of the stage "Second_Run" and, inside that
                 stage, the max-over-ranks time of the events listed in EVENTS;
  run.log      : the KSP iteration counts printed by -firedrake_1_ksp_converged_reason
                 (the last Newton solve of the run belongs to "Second_Run").
The minimum over the repeats is reported (the least disturbed run).
'''
import os
import re
import sys
import glob
import csv
from collections import defaultdict

EVENTS = ['SNESSolve', 'SNESJacobianEval', 'SNESFunctionEval', 'KSPSolve', 'PCSetUp', 'PCApply',
          'MatLUFactorNum', 'MatLUFactorSym', 'MatSolve', 'MatMult', 'KSPGMRESOrthog', 'MatCreateSubMats']


def parse_log_view(path, stage='Second_Run'):
    text = open(path).read()
    out = {}
    # ! Summary of Stages:  " 2:      Second_Run: 1.2345e+00  12.3%  ..."
    m = re.search(r'^\s*\d+:\s*' + re.escape(stage) + r':\s*([0-9.eE+-]+)', text, re.M)
    out['stage_time'] = float(m.group(1)) if m else float('nan')
    # ! Event table of the stage: from "--- Event Stage N: <stage>" to the next "--- Event Stage".
    m = re.search(r'--- Event Stage \d+: ' + re.escape(stage) + r'\s*\n(.*?)(?=\n--- Event Stage|\n-{20,}\n\n|\Z)', text, re.S)
    block = m.group(1) if m else ''
    for ev in EVENTS:
        mm = re.search(r'^' + re.escape(ev) + r'\s+(\d+)\s+[0-9.]+\s+([0-9.eE+-]+)', block, re.M)
        out[ev + '_count'] = int(mm.group(1)) if mm else 0
        out[ev] = float(mm.group(2)) if mm else 0.0
    return out


def parse_run_log(path):
    '''KSP iterations of each linear solve, grouped by Newton solve (SNES converged reason ends a group).'''
    groups, cur = [], []
    for line in open(path, errors='ignore'):
        m = re.search(r'Linear firedrake_1_ solve (converged|did not converge) due to (\w+) iterations (\d+)', line)
        if m:
            cur.append(int(m.group(3)))
        if 'Nonlinear firedrake_1_ solve' in line:
            groups.append(cur)
            cur = []
    return groups


def main(root):
    rows = []
    for d in sorted(glob.glob(os.path.join(root, '*_np*_*_rep*'))):
        m = re.match(r'(.+)_np(\d+)_(star|vanka)_rep(\d+)$', os.path.basename(d))
        if not m or not os.path.exists(os.path.join(d, 'finished')):
            continue
        name, nr, patch, rep = m.group(1), int(m.group(2)), m.group(3), int(m.group(4))
        t = parse_log_view(os.path.join(d, 'log_view.txt'))
        its = parse_run_log(os.path.join(d, 'run.log'))
        last = its[-1] if its else []
        rows.append(dict(config=name, ranks=nr, patch=patch, rep=rep, newton_steps=len(last), ksp_total=sum(last),
                         ksp_per_step='/'.join(map(str, last)), **t))
    if not rows:
        sys.exit('no finished runs under ' + root)
    keys = list(rows[0].keys())
    with open(os.path.join(root, 'timing_all_runs.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)

    best = defaultdict(dict)
    for r in rows:
        k = (r['config'], r['ranks'])
        if r['patch'] not in best[k] or r['stage_time'] < best[k][r['patch']]['stage_time']:
            best[k][r['patch']] = r
    lines = ['| config | ranks | patch | KSP its per Newton step | total KSP | Second_Run (s) | PCSetUp (s) | LU factor (s) | PCApply (s) | MatSolve (s) | time / KSP it (s) |',
             '|---|---|---|---|---|---|---|---|---|---|---|']
    ratio = ['| config | ranks | KSP its star -> vanka | time star (s) | time vanka (s) | vanka / star |', '|---|---|---|---|---|---|']
    for (cfg, nr), d in best.items():
        for patch in ('star', 'vanka'):
            if patch not in d:
                continue
            r = d[patch]
            per_it = r['stage_time'] / r['ksp_total'] if r['ksp_total'] else float('nan')
            lines.append(f"| {cfg} | {nr} | {patch} | {r['ksp_per_step']} | {r['ksp_total']} | {r['stage_time']:.2f} | {r['PCSetUp']:.2f} | "
                         f"{r['MatLUFactorNum']:.2f} | {r['PCApply']:.2f} | {r['MatSolve']:.2f} | {per_it:.3f} |")
        if 'star' in d and 'vanka' in d:
            s, v = d['star'], d['vanka']
            ratio.append(f"| {cfg} | {nr} | {s['ksp_total']} -> {v['ksp_total']} | {s['stage_time']:.2f} | {v['stage_time']:.2f} | {v['stage_time'] / s['stage_time']:.2f} |")
    n_rep = max(r['rep'] for r in rows)
    text = (f"# Timing: star vs single-column (vanka) patch\n\nMinimum over {n_rep} repeat(s) of the PETSc log stage `Second_Run` "
            "(one time step, everything compiled and set up once before).  Event times are the max over ranks, inside that stage; "
            "PCSetUp and PCApply are nested (outer PC includes the patch work).\n\n"
            + '\n'.join(ratio) + '\n\n' + '\n'.join(lines) + '\n')
    open(os.path.join(root, 'timing_summary.md'), 'w').write(text)
    print(text)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs'))
