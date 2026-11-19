'''
2026-11-18: the paper's nonlinear sweeps in the paper's shift convention, both patches, to check the paper's nonlinear
numbers (made on the Mac, data not on this machine) and its single-column comparison (supplement sm:single-column).

The nonlinear script takes the stage shift delta' of the implicit midpoint rule, which acts as twice the shift of
A_delta; the paper (since 2026-11-16, claude_small_patch/make_paper_figures_gmres6.py) therefore runs --shift C1/(2 dt)
so that delta = C1/dt as in the linear runs. The 2026-11-16/18 sweeps of this folder passed --shift C1/dt (K doubled).
This script copies every command of the back-substitution sweeps (make_backsub_jobs.py, step 1) with the shift halved;
the `# C1=` comment then holds the paper's C1. Everything else (HEAD script with the back-substitution, GMRES(6),
FGMRES right, EW, MUMPS coarse LU, 2 ranks) is unchanged.

    python make_paperconv_jobs.py      # then: ./run_pool.sh logs/joblist_paperconv.txt 12 14400
-> data/patch_gmres6/paperconv/data_nonlinear_<t>_<p>/jobs/job_<md5>/cmd.sh
'''
import glob
import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
YR2 = os.path.abspath(os.path.join(HERE, '..', '..'))
SRC = os.path.join(YR2, 'data', 'patch_gmres6')
DST = os.path.join(SRC, 'paperconv')
jobs = []
for p in ('star', 'vanka'):
    for t in ('dt', 'ar', 'dx', 'dz'):
        for cf in sorted(glob.glob(os.path.join(SRC, f'data_nonlinear_{t}_deg2_ew_gmres6_backsub_{p}', 'jobs', '*', 'cmd.sh'))):
            cmd = open(cf).read().strip()
            base, comment = cmd.split('  # ')
            s = float(re.search(r'--shift (\S+)', base).group(1))
            dt = float(re.search(r'--dt (\S+)', base).group(1))
            c1 = float(comment.split('=')[1])
            assert abs(s * dt / c1 - 1) < 1e-6, cf
            new = re.sub(r'--shift \S+', f'--shift {s / 2:.6g}', base) + f'  # C1={c1:g} (paper convention: delta = 2 shift = C1/dt)'
            j = os.path.join(DST, f'data_nonlinear_{t}_{p}', 'jobs', 'job_' + hashlib.md5(new.encode()).hexdigest()[:12])
            os.makedirs(j, exist_ok=True)
            open(os.path.join(j, 'cmd.sh'), 'w').write(new + '\n')
            jobs.append(j)
print(len(jobs), 'jobs')
dirs = sorted({os.path.dirname(os.path.dirname(j)) for j in jobs})
order = subprocess.run([sys.executable, os.path.join(HERE, 'order_jobs.py')] + dirs, capture_output=True, text=True, check=True).stdout
open(os.path.join(HERE, 'logs', 'joblist_paperconv.txt'), 'w').write(order)
print(len(order.split()), 'jobs in logs/joblist_paperconv.txt')
