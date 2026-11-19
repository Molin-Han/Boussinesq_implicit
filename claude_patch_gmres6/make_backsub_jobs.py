'''
Follow-up of 2026-11-18 (two steps, both with the nonlinear code at HEAD, whose auxiliary-operator preconditioner applies
the pressure back-substitution by default since commit 6c4b0cf, 2026-11-17; --no_backsub gives the preconditioner of the
2026-11-16 sweeps):

  step 1  the paper's four nonlinear sweeps for both patches, with exactly the commands of the 2026-11-16 sweeps
          (Bash_Generator/bash_paper_patch.py, MUMPS coarse LU, 2 ranks)
          -> data/patch_gmres6/data_nonlinear_<t>_deg2_ew_gmres6_backsub_<p>/   (made by run_claude_paper.sh -n)
  step 2  every point of those sweeps with K = delta/delta_thr < 1, both patches, with SuperLU_DIST on the multigrid
          coarse level (deterministic, see figures/60_mpi/README.md), run twice
          -> data/patch_gmres6/backsub_superlu/data_nonlinear_<t>_<p>_r<R>/

    python make_backsub_jobs.py        # then: ./run_pool.sh logs/joblist_backsub.txt 12 14400

The pilot that checked the set-up (HEAD back-substitution = the 2026-10-16 --pc_fix script; --no_backsub = the old counts;
two SuperLU_DIST repeats bitwise identical) is in data/patch_gmres6/backsub/pilot/.
'''
import glob
import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
YR2 = os.path.abspath(os.path.join(REPO, '..'))
NEW = os.path.join(YR2, 'data', 'patch_gmres6')
sys.path.insert(0, REPO)
import claude_paper_shift as cps   # noqa: E402

LX, N_FREQ = 3.0e5, 0.01
SLU = '-firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_pc_factor_mat_solver_type superlu_dist'
TESTS = ('dt', 'ar', 'dx', 'dz')
PATCHES = ('star', 'vanka')

# step 1: job folders through the standard driver (also writes sweeps/patch_gmres6/test_paper_*_backsub_<p>.sh)
for p in PATCHES:
    subprocess.run(['./run_claude_paper.sh', '-n', '-g', 'Bash_Generator/bash_paper_patch.py', '-D', 'data', '-S', 'test_paper',
                    '-C', 'patch_gmres6', '-m', 'nonlinear', '-t', ' '.join(TESTS), '-s', f'_deg2_ew_gmres6_backsub_{p}',
                    '-x', f'--patch {p}'], cwd=REPO, check=True, stdout=subprocess.DEVNULL)

jobs1, jobs2 = [], []
for p in PATCHES:
    for t in TESTS:
        d = os.path.join(NEW, f'data_nonlinear_{t}_deg2_ew_gmres6_backsub_{p}')
        for cf in sorted(glob.glob(os.path.join(d, 'jobs', '*', 'cmd.sh'))):
            jobs1.append(os.path.dirname(cf))
            cmd = open(cf).read().strip()
            g = lambda k: float(re.search(rf'--{k} (\S+)', cmd).group(1))   # noqa: E731
            K = g('shift') / cps.shift_threshold(g('dt'), LX, g('height'), N_FREQ)
            if K >= 1:
                continue
            base, comment = cmd.split('  # ')
            new = f'{base} {SLU}  # {comment}'
            for r in (1, 2):
                j = os.path.join(NEW, 'backsub_superlu', f'data_nonlinear_{t}_{p}_r{r}', 'jobs',
                                 'job_' + hashlib.md5(new.encode()).hexdigest()[:12])
                os.makedirs(j, exist_ok=True)
                open(os.path.join(j, 'cmd.sh'), 'w').write(new + '\n')
                jobs2.append(j)
assert '--no_backsub' not in ''.join(open(os.path.join(j, 'cmd.sh')).read() for j in jobs1 + jobs2)
print(len(jobs1), 'step-1 jobs,', len(jobs2), 'step-2 jobs')
dirs = sorted({os.path.dirname(os.path.dirname(j)) for j in jobs1 + jobs2})
order = subprocess.run([sys.executable, os.path.join(HERE, 'order_jobs.py')] + dirs, capture_output=True, text=True, check=True).stdout
open(os.path.join(HERE, 'logs', 'joblist_backsub.txt'), 'w').write(order)
print(len(order.split()), 'jobs in logs/joblist_backsub.txt')
