'''
Noise check for the MPI study: is the run-to-run variation of the nonlinear single-column counts at small K caused by the
parallel MUMPS coarse LU?  Two nonlinear single-column points of make_mpi_jobs.py (P3: dt 5e4 s, C1 1e-5; P9: Nz 80,
dt 5000 s), run
  * serial (np = 1, PETSc's own LU on the coarse level) twice more  -> data/patch_gmres6/mpi/noise_check/np01_petscLU_r{2,3}
    (r1 is the np = 1 run of the MPI study, data/patch_gmres6/mpi/nonlinear_vanka/np01_r1)
  * on 2 ranks with SuperLU_DIST instead of MUMPS on the coarse level, three times
                                                                    -> data/patch_gmres6/mpi/noise_check/np02_superlu_r{1,2,3}
and compared (counts and bitwise output files) with the MUMPS runs at np = 2 (np02_r{1,2,3} of the MPI study).

    python make_noise_jobs.py      # then: ./run_pool.sh logs/joblist_noise.txt 4
'''
import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.join(REPO, '..', 'data', 'patch_gmres6', 'mpi', 'noise_check')
GEN = os.path.join(REPO, 'Bash_Generator', 'bash_paper_patch.py')
POINTS = [('P3', ['--test', 'dt', '--C1_list', '1e-05', '--dts', '50000.0']),
          ('P9', ['--test', 'dz', '--nzs', '80', '--dts', '5000.0'])]
SUPERLU = '-firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_pc_factor_mat_solver_type superlu_dist'
jobs = []
for variant, np_, extra, reps in (('np01_petscLU', 1, None, (2, 3)), ('np02_superlu', 2, SUPERLU, (1, 2, 3))):
    for pid, args in POINTS:
        a = [sys.executable, '-I', '-B', GEN, '--mode', 'nonlinear', '--patch', 'vanka', '--ncpus', str(np_)] + args
        if extra:
            a.append('--extra=' + extra)
        cmds = [l for l in subprocess.run(a, capture_output=True, text=True, check=True).stdout.splitlines() if l.startswith('mpiexec')]
        assert len(cmds) == 1
        cmd = re.sub(r'(\.\./)+Boussinesq_implicit/', REPO + '/', cmds[0])
        cmd = re.sub(r'^mpiexec -n', 'mpiexec --bind-to none -n', cmd)
        for r in reps:
            j = os.path.join(os.path.abspath(ROOT), f'{variant}_r{r}', 'jobs', 'job_' + hashlib.md5(cmd.encode()).hexdigest()[:12])
            os.makedirs(j, exist_ok=True)
            open(os.path.join(j, 'cmd.sh'), 'w').write(cmd + '\n')
            jobs.append(j)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
open(os.path.join(HERE, 'logs', 'joblist_noise.txt'), 'w').write('\n'.join(jobs) + '\n')
print(len(jobs), 'jobs;', jobs[0])
