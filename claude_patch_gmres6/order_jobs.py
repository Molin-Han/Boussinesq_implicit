'''
Order job directories longest-first so that a single xargs pool finishes without a long tail.

    python order_jobs.py <data dir> [<data dir> ...] > joblist.txt

The expected run time of a job is looked up from the 2026-10-16 runs of the same parameters
(data/paper_robustness/data_<mode>_<test>_deg2[_ew]_<patch>/jobs/*/status, 2 ranks); jobs without a
match get a large default so that they start early.  The rank count (mpiexec -n) scales the estimate.
'''
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
YR2 = os.path.abspath(os.path.join(HERE, '..', '..'))
KEYS = ('dt', 'nx', 'nz', 'height', 'shift')


def params(cmd):
    p = {k: float(m.group(1)) for k in KEYS if (m := re.search(rf'--{k} (\S+)', cmd))}
    p['code'] = 'nonlinear' if 'Nonlinear_Boussinesq_Irk_SC' in cmd else 'linear'
    p['patch'] = 'vanka' if re.search(r'--patch[ =]vanka', cmd) else 'star'
    m = re.search(r'-n (\d+) ', cmd)
    p['np'] = int(m.group(1)) if m else 2
    return p


def key(p):
    return (p['code'], p['patch'], p['dt'], p['nx'], p['nz'], p['height'], round(p['shift'], 18))


known = {}
for st in glob.glob(f'{YR2}/data/paper_robustness/data_*_deg2*_*/jobs/*/status'):
    try:
        secs = int(open(st).read().split()[1])
        p = params(open(os.path.join(os.path.dirname(st), 'cmd.sh')).read())
    except (OSError, IndexError, ValueError, AttributeError):
        continue
    known[key(p)] = max(secs, known.get(key(p), 0))

jobs = []
for d in sys.argv[1:]:
    for j in sorted(glob.glob(os.path.join(os.path.abspath(d), 'jobs', 'job_*'))):
        p = params(open(os.path.join(j, 'cmd.sh')).read())
        est = known.get(key(p), 3000) * (2.0 / p['np']) ** 0.5
        jobs.append((est, j))
for est, j in sorted(jobs, reverse=True):
    print(j)
print(f'{len(jobs)} jobs, estimated {sum(e for e, _ in jobs) / 3600:.1f} job-hours (2-rank-equivalent)', file=sys.stderr)
