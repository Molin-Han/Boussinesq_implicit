'''
Job folders for the MPI-robustness study: the same solver configurations as the main GMRES(6) sweeps
(Bash_Generator/bash_paper_patch.py, delta = C1/dt), at selected parameter points, run on np = 1, 2, 4, 8, 16
MPI ranks, with repeats at np = 2 and np = 8 to measure the run-to-run noise (parallel MUMPS coarse LU).

    python make_mpi_jobs.py            # writes ../../data/patch_gmres6/mpi/<code>_<patch>/np<NN>_r<R>/jobs/job_<md5>/cmd.sh
                                       # and points.csv, and prints one joblist per np to logs/joblist_mpi_np<NN>.txt

Each command is exactly what bash_paper_patch.py prints for that point with --ncpus <np>, made absolute and
given --bind-to none the same way as run_claude_paper.sh does; run them with run_pool.sh.
'''
import csv
import hashlib
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..'))
YR2 = os.path.abspath(os.path.join(REPO, '..'))
ROOT = os.path.join(YR2, 'data', 'patch_gmres6', 'mpi')
GEN = os.path.join(REPO, 'Bash_Generator', 'bash_paper_patch.py')
PY = sys.executable

NPS = [1, 2, 4, 8, 16]
REPEATS = {2: 3, 8: 3}          # np -> number of runs (base + repeats); 1 run otherwise

# point id -> (description, generator arguments) per code.  dt = 50 s unless the point is about dt.
POINTS = {
    'linear': [
        ('P1', 'nominal (20x20 coarse, Lz 4 km, dt 50 s)', ['--test', 'dx', '--nxs', '20', '--dts', '50.0']),
        ('P2', 'dt 1 s (C1 1e-3)', ['--test', 'dt', '--C1_list', '0.001', '--dts', '1.0']),
        ('P3', 'dt 1e4 s (C1 1e-3, smallest shift)', ['--test', 'dt', '--C1_list', '0.001', '--dts', '10000.0']),
        ('P4', 'Lz 12 km', ['--test', 'ar', '--heights', '12000', '--dts', '50.0']),
        ('P5', 'Lz 0.25 km', ['--test', 'ar', '--heights', '250', '--dts', '50.0']),
        ('P6', 'Nx 80', ['--test', 'dx', '--nxs', '80', '--dts', '50.0']),
        ('P7', 'Nx 10', ['--test', 'dx', '--nxs', '10', '--dts', '50.0']),
        ('P8', 'Nz 80', ['--test', 'dz', '--nzs', '80', '--dts', '50.0']),
        ('P9', 'stress: Lz 0.25 km, dt 1000 s', ['--test', 'ar', '--heights', '250', '--dts', '1000.0']),
    ],
    'nonlinear': [
        ('P1', 'nominal (20x20 coarse, Lz 4 km, dt 50 s)', ['--test', 'dx', '--nxs', '20', '--dts', '50.0']),
        ('P2', 'dt 1 s (C1 1e-5)', ['--test', 'dt', '--C1_list', '1e-05', '--dts', '1.0']),
        ('P3', 'dt 5e4 s (C1 1e-5, smallest shift)', ['--test', 'dt', '--C1_list', '1e-05', '--dts', '50000.0']),
        ('P4', 'Lz 12 km', ['--test', 'ar', '--heights', '12000', '--dts', '50.0']),
        ('P5', 'Lz 1 km', ['--test', 'ar', '--heights', '1000', '--dts', '50.0']),
        ('P6', 'Nx 80', ['--test', 'dx', '--nxs', '80', '--dts', '50.0']),
        ('P7', 'Nx 10', ['--test', 'dx', '--nxs', '10', '--dts', '50.0']),
        ('P8', 'Nz 80', ['--test', 'dz', '--nzs', '80', '--dts', '50.0']),
        ('P9', 'stress: Nz 80, dt 5000 s', ['--test', 'dz', '--nzs', '80', '--dts', '5000.0']),
    ],
}


def commands(code, patch, np_, args):
    extra = ['--patch', patch, '--ncpus', str(np_)]
    if code == 'linear':
        extra.append('--extra=--mg_smoother=gmres')
    out = subprocess.run([PY, '-I', '-B', GEN, '--mode', code] + args + extra, capture_output=True, text=True, check=True).stdout
    return [l for l in out.splitlines() if l.startswith('mpiexec')]


rows, lists = [], {n: [] for n in NPS}
for code, pts in POINTS.items():
    for pid, desc, args in pts:
        nx = int(args[args.index('--nxs') + 1]) if '--nxs' in args else 20
        for patch in ('star', 'vanka'):
            for np_ in NPS:
                if np_ > nx:
                    continue        # the coarse mesh needs at least one cell per rank
                cmds = commands(code, patch, np_, args)
                assert len(cmds) == 1, (code, pid, cmds)
                cmd = re.sub(r'(\.\./)+Boussinesq_implicit/', REPO + '/', cmds[0])
                cmd = re.sub(r'^mpiexec -n', 'mpiexec --bind-to none -n', cmd)
                for r in range(1, REPEATS.get(np_, 1) + 1):
                    d = os.path.join(ROOT, f'{code}_{patch}', f'np{np_:02d}_r{r}')
                    h = hashlib.md5(cmd.encode()).hexdigest()[:12]
                    j = os.path.join(d, 'jobs', f'job_{h}')
                    os.makedirs(j, exist_ok=True)
                    with open(os.path.join(j, 'cmd.sh'), 'w') as f:
                        f.write(cmd + '\n')
                    lists[np_].append(j)
                    rows.append(dict(code=code, point=pid, description=desc, patch=patch, np=np_, repeat=r,
                                     folder=os.path.relpath(d, ROOT), job=os.path.basename(j), command=cmd))
os.makedirs(ROOT, exist_ok=True)
with open(os.path.join(ROOT, 'points.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
for np_, js in lists.items():
    with open(os.path.join(HERE, 'logs', f'joblist_mpi_np{np_:02d}.txt'), 'w') as f:
        f.write('\n'.join(js) + '\n')
    print(f'np={np_:2d}: {len(js)} jobs')
print(f'total {len(rows)} jobs; table {os.path.join(ROOT, "points.csv")}')
