"""
Item 4: iteration counts and residual histories of the solver script for different numbers
of MPI ranks (and for the committed version of the script).

Reads the run.log files written by run_case.sh (runs/<label>/run.log) and the two logs
that were in the folder at the start (run_star_nominal_2ranks.log, run_vanka_nominal_2ranks.log).
For every log the KSP solves are listed in the order in which they happen:
    time step 1:            Newton 1, Newton 2
    time step 2 (warm-up):  Newton 1, Newton 2
    time step 2 (2nd run):  Newton 1, Newton 2     <- this is what snes_ksp_cum_*.out contains
and the residual histories of two logs are compared iteration by iteration.

Usage: python 04_compare_residual_histories.py            (writes 04_residual_histories.txt)
"""
import os
import re
import sys
import glob

HERE = os.path.dirname(os.path.abspath(__file__))

KSP_RE = re.compile(r"^\s*(\d+) KSP unpreconditioned resid norm ([0-9.eE+-]+) true resid norm ([0-9.eE+-]+)")
SNES_RE = re.compile(r"^\s*(\d+) SNES Function norm ([0-9.eE+-]+)")
REASON_RE = re.compile(r"Linear firedrake_1_ solve converged due to (\w+) iterations (\d+)")
REASON_BAD_RE = re.compile(r"Linear firedrake_1_ solve did not converge due to (\w+) iterations (\d+)")


def parse(path):
    """Return a list of KSP solves: dict(res=[...], reason=str, snes0=float)."""
    solves = []
    cur = None
    snes = []
    for line in open(path, errors='replace'):
        m = SNES_RE.match(line)
        if m:
            snes.append((int(m.group(1)), float(m.group(2))))
            continue
        m = KSP_RE.match(line)
        if m:
            it = int(m.group(1))
            if it == 0:
                cur = dict(res=[], reason=None)
                solves.append(cur)
            cur['res'].append(float(m.group(2)))
            continue
        m = REASON_RE.search(line) or REASON_BAD_RE.search(line)
        if m and cur is not None and cur['reason'] is None:
            cur['reason'] = m.group(1)
    return solves, snes


def label_of(path):
    if os.path.basename(path) == 'run.log':
        return os.path.basename(os.path.dirname(path))
    return os.path.basename(path).replace('.log', '') + ' (provided)'


logs = sorted(glob.glob(os.path.join(HERE, 'runs', '*', 'run.log')))
logs = [p for p in logs if not os.path.basename(os.path.dirname(p)).startswith('patch_dofs')]
logs += [os.path.join(HERE, f) for f in ('run_star_nominal_2ranks.log', 'run_vanka_nominal_2ranks.log')
         if os.path.exists(os.path.join(HERE, f))]
data = {}
for p in logs:
    s, n = parse(p)
    if s:
        data[label_of(p)] = (s, n)

lines = []
out = lines.append
NAMES = ['step1 Newton1', 'step1 Newton2', 'step2 warm-up Newton1', 'step2 warm-up Newton2',
         'step2 2nd-run Newton1', 'step2 2nd-run Newton2']

out("FGMRES iterations of every linear solve (iterations, convergence reason)")
out(f"{'run':38s} " + " ".join(f"{n:>22s}" for n in NAMES) + "   cumulative 2nd run")
for lab in sorted(data):
    s, _ = data[lab]
    cells = []
    for k in range(len(NAMES)):
        if k < len(s):
            cells.append(f"{len(s[k]['res']) - 1:3d} {str(s[k]['reason'] or '?').replace('CONVERGED_', ''):>18s}")
        else:
            cells.append(f"{'-':>22s}")
    cum = ""
    if len(s) >= 6:
        a = len(s[4]['res']) - 1
        b = len(s[5]['res']) - 1
        cum = f"0/{a}/{a + b}"
    out(f"{lab:38s} " + " ".join(cells) + f"   {cum}")


def compare(a, b, which):
    if a not in data or b not in data:
        return
    sa, sb = data[a][0], data[b][0]
    out("")
    out(f"relative difference of the residual norms, {a}  versus  {b}")
    for k in which:
        if k >= len(sa) or k >= len(sb):
            continue
        ra, rb = sa[k]['res'], sb[k]['res']
        n = min(len(ra), len(rb))
        out(f"  {NAMES[k]}: iterations {len(ra) - 1} / {len(rb) - 1}")
        out("     it      " + f"{a[:24]:>24s} {b[:24]:>24s}   rel.diff   |diff|/r0")
        for i in range(n):
            out(f"     {i:3d}     {ra[i]:24.12e} {rb[i]:24.12e}   {abs(ra[i] - rb[i]) / abs(ra[i]):8.1e}   "
                f"{abs(ra[i] - rb[i]) / ra[0]:8.1e}")


pairs = [('star_np1', 'committedHEAD_np1'), ('star_np2', 'committedHEAD_np2'),
         ('star_np2', 'run_star_nominal_2ranks (provided)'), ('vanka_np2', 'run_vanka_nominal_2ranks (provided)'),
         ('star_np1', 'star_np2'), ('vanka_np1', 'vanka_np2'),
         ('star_np1', 'star_np4'), ('vanka_np1', 'vanka_np4'),
         ('star_np2', 'star_np4'), ('vanka_np2', 'vanka_np4'),
         ('vanka_np1', 'vanka_np1_repeat'), ('vanka_np2', 'vanka_np2_repeat')]
for a, b in pairs:
    compare(a, b, [0, 4, 5])

text = "\n".join(lines)
print(text)
with open(os.path.join(HERE, '04_residual_histories.txt'), 'w') as f:
    f.write(text + "\n")
