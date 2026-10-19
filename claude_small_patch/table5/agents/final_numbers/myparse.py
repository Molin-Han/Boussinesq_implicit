"""Independent parser of runs/*: outer counts, Official-Run PCApply/KSPSolve, patch sizes, nprocs, status, cmd."""
import os, re, glob, json, sys
RUNS = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/runs'

def stages(logfile):
    out = {}
    if not os.path.exists(logfile):
        return out
    cur = None
    with open(logfile) as f:
        for line in f:
            m = re.match(r'--- Event Stage (\d+): (.*)$', line.strip())
            if m:
                cur = m.group(2).strip()
                out.setdefault(cur, {})
                continue
            if cur is not None:
                p = line.split()
                if len(p) > 4 and p[0] in ('PCApply', 'KSPSolve', 'PCSetUp', 'SNESSolve', 'MatMult') and p[0] not in out[cur]:
                    try:
                        out[cur][p[0]] = (int(p[1]), float(p[3]))
                    except ValueError:
                        pass
    return out

def nprocs(logfile):
    if not os.path.exists(logfile):
        return None
    with open(logfile) as f:
        for line in f:
            m = re.search(r'with (\d+) process', line)
            if m:
                return int(m.group(1))
    return None

def outer(runlog):
    its = []
    if not os.path.exists(runlog):
        return its
    with open(runlog, errors='replace') as f:
        for line in f:
            m = re.search(r'Linear firedrake_2_ solve (converged|did not converge) due to (\S+) iterations (\d+)', line)
            if m:
                its.append((int(m.group(3)), m.group(2)))
    return its

def patches(f):
    s = set()
    if not os.path.exists(f):
        return s
    with open(f, errors='replace') as fh:
        for line in fh:
            m = re.search(r'Minimum / average / maximum patch sizes : (\d+) / ([\d.]+) / (\d+)', line)
            if m:
                s.add(int(m.group(3)))
    return s

rows = {}
for d in sorted(glob.glob(os.path.join(RUNS, '*'))):
    tag = os.path.basename(d)
    cmd = open(os.path.join(d, 'cmd.sh')).read().strip() if os.path.exists(os.path.join(d, 'cmd.sh')) else ''
    mnp = re.search(r'-n (\d+)', cmd)
    st = open(os.path.join(d, 'status')).read().split() if os.path.exists(os.path.join(d, 'status')) else []
    sg = stages(os.path.join(d, 'log.txt'))
    off = sg.get('Official-Run', {})
    wu = sg.get('Warm-up-solve', {})
    its = outer(os.path.join(d, 'run.log'))
    pt = patches(os.path.join(d, 'run.log')) | patches(os.path.join(d, 'slice3D.txt'))
    rows[tag] = dict(cmd=cmd, np_cmd=int(mnp.group(1)) if mnp else None, np_log=nprocs(os.path.join(d, 'log.txt')),
                     status=st, its=[i for i, r in its], reasons=[r for i, r in its],
                     off_pc=off.get('PCApply'), off_ksp=off.get('KSPSolve'), wu_pc=wu.get('PCApply'), wu_ksp=wu.get('KSPSolve'),
                     patch=sorted(pt), files=sorted(os.listdir(d)))
json.dump(rows, open('/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/table5/agents/final_numbers/myparse.json', 'w'), indent=1)
for t, r in rows.items():
    print(t, 'np', r['np_cmd'], r['np_log'], 'st', r['status'], 'its', '/'.join(map(str, r['its'])), set(r['reasons']),
          'offPC', r['off_pc'], 'offKSP', r['off_ksp'], 'patch', r['patch'])
