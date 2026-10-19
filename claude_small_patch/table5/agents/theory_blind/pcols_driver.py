"""Save, for every pressure dof of LB_time_slice.py's finest space, the x-centre of its column
(cell) and its own (x,z): runs the script's set-up only (no solve). Same argv as the dense run."""
import sys
import numpy as np
REPO = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit'
sys.path.insert(0, REPO)
src = open(REPO + '/LB_time_slice.py').read()
head = src.split('# Time Stepping')[0]
sys.argv = ['LB_time_slice.py'] + sys.argv[1:]
exec(compile(head, REPO + '/LB_time_slice.py', 'exec'), globals())
DG0 = FunctionSpace(mesh, 'DG', 0)
xc = Function(Pressure).interpolate(Function(DG0).interpolate(x)).dat.data_ro.copy()
zc = Function(Pressure).interpolate(Function(DG0).interpolate(z)).dat.data_ro.copy()
xp = Function(Pressure).interpolate(x).dat.data_ro.copy()
zp = Function(Pressure).interpolate(z).dat.data_ro.copy()
np.savez('pcols.npz', xc=xc, zc=zc, xp=xp, zp=zp)
print('saved pcols.npz', xc.shape)
