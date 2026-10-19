"""Dense diagnostics of the outer preconditioner of LB_time_slice.py (serial, tiny meshes only).

Executes the UNMODIFIED LB_time_slice.py up to its time loop (same argv, same options, same
firedrake_2_ prefix), registers a KSP monitor that -- at outer iteration 0, i.e. inside
snes.solve where all DM hooks are active -- applies the outer PC to every unit vector, and stores
  Pinv (dense P^{-1}), A (outer Amat = unshifted J), Ad (outer Pmat = shifted Jp), field offsets,
  and the (x,z) coordinates of the pressure dofs
in <cwd>/dense.npz.  Usage: python spectra_driver.py <LB_time_slice.py args>  (with --tmax 100).
"""
import sys
import numpy as np
REPO = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit'
sys.path.insert(0, REPO)
src = open(REPO + '/LB_time_slice.py').read()
head = src.split('# Time Stepping')[0]
sys.argv = ['LB_time_slice.py'] + sys.argv[1:]
exec(compile(head, REPO + '/LB_time_slice.py', 'exec'), globals())

from scipy.sparse import csr_matrix
store = {}

def dense(M):
    ai, aj, av = M.getValuesCSR()
    return csr_matrix((av, aj, ai), shape=M.getSize()).toarray()

def mon(ksp, it, rnorm):
    if it != 0 or 'Pinv' in store:
        return
    pc = ksp.getPC()
    Amat, Pmat = ksp.getOperators()
    n = Amat.getSize()[0]
    x, y = Amat.createVecs()
    Pinv = np.empty((n, n))
    for j in range(n):
        x.zeroEntries(); x.setValue(j, 1.0); x.assemble()
        pc.apply(x, y)
        Pinv[:, j] = y.getArray(readonly=True)
    store['Pinv'] = Pinv
    store['A'] = dense(Amat)
    store['Ad'] = dense(Pmat)

nsolver.snes.ksp.setMonitor(mon)
Unp1.assign(Un)
nsolver.solve()
print('outer iterations', nsolver.snes.ksp.getIterationNumber())
sizes = [Wi.dim() for Wi in W]
xp = Function(Pressure).interpolate(x).dat.data_ro.copy()
zp = Function(Pressure).interpolate(z).dat.data_ro.copy()
np.savez('dense.npz', Pinv=store['Pinv'], A=store['A'], Ad=store['Ad'], sizes=np.array(sizes),
         xp=xp, zp=zp, shift=args.shift, dt=args.dt)
print('saved dense.npz, sizes', sizes)
