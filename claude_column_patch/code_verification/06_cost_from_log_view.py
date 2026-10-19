"""
Item 5 (fairness): iteration counts are not cost.  This script reads the PETSc -log_view files
written by run_case.sh (runs/<label>/log_view.txt) and extracts, for the stage 'Second_Run'
(the measured second time step), the quantities that tell what one outer FGMRES iteration costs
with the star and with the vanka patch.

CAUTION.  The wall-clock times were measured at nice 19 on a machine that was fully loaded by
a parameter sweep, so they are indicative only.  The flop counts (PETSc operations only: the
matrix-free Jacobian and the assembly are PyOP2 kernels and only partly counted) and the call
counts are deterministic.

Usage: python 06_cost_from_log_view.py       -> 06_cost_from_log_view.txt
"""
import os
import re
import glob

HERE = os.path.dirname(os.path.abspath(__file__))
EVENTS = ["SNESSolve", "KSPSolve", "PCSetUp", "PCApply", "PCApplyOnBlocks", "MatSolve", "MatLUFactorNum",
          "MatLUFactorSym", "MatCreateSubMats", "MatMult", "KSPGMRESOrthog", "SNESFunctionEval",
          "SNESJacobianEval", "ParLoopExecute"]


def parse(path):
    txt = open(path, errors='replace').read()
    m = re.search(r"--- Event Stage \d+: Second_Run(.*?)(--- Event Stage|\n-{60,}\n)", txt, re.S)
    if not m:
        return None
    block = m.group(1)
    ev = {}
    for line in block.splitlines():
        parts = line.split()
        if len(parts) > 6 and parts[0] in EVENTS:
            try:
                ev[parts[0]] = dict(count=int(parts[1]), time=float(parts[3]), flops=float(parts[5]))
            except ValueError:
                pass
    return ev


def its_of(d):
    f = glob.glob(os.path.join(d, "snes_ksp_cum_*.out"))
    if not f:
        return None
    v = [int(float(x)) for x in open(f[0]).read().split()]
    return v


rows = []
for d in sorted(glob.glob(os.path.join(HERE, "runs", "*"))):
    lv = os.path.join(d, "log_view.txt")
    if not os.path.exists(lv):
        continue
    ev = parse(lv)
    its = its_of(d)
    if ev is None or its is None:
        continue
    rows.append((os.path.basename(d), its, ev, open(os.path.join(d, "walltime.txt")).read().strip()))

lines = []
out = lines.append
out("Stage 'Second_Run' (second time step, all Newton steps).  time = max over ranks (s), flops = max over ranks.")
out("")
hdr = f"{'run':24s} {'cum. its':>10s} {'outer its':>9s}"
for e in ["SNESSolve", "PCSetUp", "PCApply", "MatSolve", "MatLUFactorNum", "MatMult"]:
    hdr += f" {e + ' t':>17s} {e + ' fl':>17s}"
out(hdr)
for lab, its, ev, wt in rows:
    n = its[-1]
    s = f"{lab:24s} {'/'.join(map(str, its)):>10s} {n:9d}"
    for e in ["SNESSolve", "PCSetUp", "PCApply", "MatSolve", "MatLUFactorNum", "MatMult"]:
        if e in ev:
            s += f" {ev[e]['time']:17.3f} {ev[e]['flops']:17.3e}"
        else:
            s += f" {'-':>17s} {'-':>17s}"
    out(s)
out("")
out("per outer FGMRES iteration (PCApply = one application of the whole preconditioner):")
out(f"{'run':24s} {'PCApply calls':>13s} {'t(PCApply)/call':>16s} {'flops(PCApply)/call':>20s} {'patch solves/call':>18s} "
    f"{'t(SNESSolve)/it':>16s} {'t(PCSetUp)':>11s} {'wall clock of the whole run':>30s}")
for lab, its, ev, wt in rows:
    n = its[-1]
    if 'PCApply' not in ev:
        continue
    pa = ev['PCApply']
    out(f"{lab:24s} {pa['count']:13d} {pa['time'] / pa['count']:16.3f} {pa['flops'] / pa['count']:20.3e} "
        f"{ev['MatSolve']['count'] / pa['count']:18.1f} {ev['SNESSolve']['time'] / max(n, 1):16.3f} "
        f"{ev['PCSetUp']['time']:11.3f}   {wt}")
text = "\n".join(lines)
print(text)
with open(os.path.join(HERE, "06_cost_from_log_view.txt"), "w") as f:
    f.write(text + "\n")
