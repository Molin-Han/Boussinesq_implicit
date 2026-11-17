'''
Diagnosis of the two runs of the MPI study that aborted: nonlinear code, P6 (Nx = 80, dt = 50 s, C1 = 1e-5) at np = 16,
both patches (data/patch_gmres6/mpi/nonlinear_{star,vanka}/np16_r1).  Both stopped in the first Newton step with
    KSPConvergedRateView -> KSPComputeConvergenceRate -> PetscLinearRegression: LAPACKgetrf failed with INFO = 2,
i.e. the outer KSP returned after 0 iterations and the ksp_converged_rate diagnostic (set in the script's option dict,
which the command line cannot override) failed on a residual history of length 1.  These reruns print WHY the KSP
stopped (converged-reason views are printed before the rate view), and test whether the run passes when the coarse LU
uses SuperLU_DIST instead of MUMPS:

  data/patch_gmres6/mpi/np16_check/<patch>_mumps/    the failing command + converged-reason views
  data/patch_gmres6/mpi/np16_check/<patch>_superlu/  the same with -..._mg_coarse_pc_factor_mat_solver_type superlu_dist
  data/patch_gmres6/mpi/np16_check/<patch>_mumps_icntl14/  MUMPS with 200 % extra working space (ICNTL(14) = 200,
                                                     MUMPS default 20 %) and its error messages on (ICNTL(4) = 1)

The first two variants showed that the MUMPS factorisation of the coarse level fails (PC failed due to FACTOR_OUTMEMORY),
hence the third.

    python make_np16_check_jobs.py      # then: ./run_pool.sh logs/joblist_np16_check.txt 2
'''
import glob
import hashlib
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MPI = os.path.abspath(os.path.join(HERE, '..', '..', 'data', 'patch_gmres6', 'mpi'))
REASON = ' '.join('-firedrake_1_' + p + 'ksp_converged_reason' for p in (
    '', 'shiftedschurpc_fieldsplit_0_', 'shiftedschurpc_fieldsplit_1_', 'shiftedschurpc_fieldsplit_1_mg_coarse_',
    'shiftedschurpc_fieldsplit_1_mg_levels_'))
SUPERLU = '-firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_pc_factor_mat_solver_type superlu_dist'
ICNTL14 = '-firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_mat_mumps_icntl_14 200 -firedrake_1_shiftedschurpc_fieldsplit_1_mg_coarse_mat_mumps_icntl_4 1'
jobs = []
for patch in ('star', 'vanka'):
    src = [f for f in glob.glob(os.path.join(MPI, f'nonlinear_{patch}', 'np16_r1', 'jobs', '*', 'cmd.sh'))
           if '--nx 80 ' in open(f).read()]
    assert len(src) == 1, src
    base, comment = open(src[0]).read().strip().split('  # ')
    for variant, extra in (('mumps', REASON), ('superlu', SUPERLU + ' ' + REASON), ('mumps_icntl14', ICNTL14 + ' ' + REASON)):
        cmd = f'{base} {extra}  # {comment}'
        j = os.path.join(MPI, 'np16_check', f'{patch}_{variant}', 'jobs', 'job_' + hashlib.md5(cmd.encode()).hexdigest()[:12])
        os.makedirs(j, exist_ok=True)
        open(os.path.join(j, 'cmd.sh'), 'w').write(cmd + '\n')
        jobs.append(j)
open(os.path.join(HERE, 'logs', 'joblist_np16_check.txt'), 'w').write('\n'.join(jobs) + '\n')
open(os.path.join(HERE, 'logs', 'joblist_np16_check_icntl14.txt'), 'w').write('\n'.join(j for j in jobs if '_icntl14' in j) + '\n')
print(len(jobs), 'jobs'); print('\n'.join(jobs))
