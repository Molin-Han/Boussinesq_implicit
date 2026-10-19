'''
Is the prolongation that Firedrake uses in the multigrid of field 1 the natural embedding
of the coarse space into the fine space (so that the rediscretised coarse operator is the
Galerkin operator, as assumed in emulator.py)?

Test: for random coarse vectors x_c, compare the energies x_c^T A_c x_c and
(P x_c)^T A_f (P x_c), and the L2 norms, for the velocity space and the buoyancy space,
with P = firedrake prolong (the TransferManager default used by the solver).

Run (serial): python verify_transfer_firedrake.py      (writes verify_transfer_firedrake.out)
'''
import sys
import os
import numpy as np
sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
from firedrake import *          # noqa
import utils

here = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(here, 'verify_transfer_firedrake.out'), 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.write(s + '\n')
    out.flush()


nx, nz, length, height, deg = 8, 6, 8 * 15000.0, 4000.0, 2
c, m = 2.279727e9, 1.0625
distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
base = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(base, 1)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * 2, extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mc, mf = new_mh[0], new_mh[1]
kz = utils.k()
rng = np.random.default_rng(1)


def energy_forms(V):
    u, w = TrialFunction(V), TestFunction(V)
    if V.value_shape == ():
        return inner(u, w) * dx, inner(u, w) * dx
    a = (inner(u, w) + Constant(m - 1.0) * inner(u, kz) * inner(w, kz) + Constant(c) * div(u) * div(w)) * dx
    return a, inner(u, w) * dx


for name, make in (('velocity V_2D', lambda msh: utils.extrude_RT(msh, k=deg)),
                   ('buoyancy W_theta', lambda msh: utils.W_theta(msh, k=deg)),
                   ('u_y DG', lambda msh: FunctionSpace(msh, 'DG', deg - 1))):
    Vc, Vf = make(mc), make(mf)
    ac, mcf = energy_forms(Vc)
    af, mff = energy_forms(Vf)
    bcs_c = bcs_f = []
    if Vc.value_shape != ():
        bcs_c = [DirichletBC(Vc, as_vector([0., 0., 0.]), s) for s in ("top", "bottom")]
        bcs_f = [DirichletBC(Vf, as_vector([0., 0., 0.]), s) for s in ("top", "bottom")]
    worst_a, worst_m = 0.0, 0.0
    for trial in range(5):
        xc = Function(Vc)
        xc.dat.data[:] = rng.standard_normal(xc.dat.data.shape)
        for bc in bcs_c:
            bc.zero(xc)
        xf = Function(Vf)
        prolong(xc, xf)
        ec = assemble(action(action(ac, xc), xc))
        ef = assemble(action(action(af, xf), xf))
        nc = assemble(action(action(mcf, xc), xc))
        nf = assemble(action(action(mff, xf), xf))
        worst_a = max(worst_a, abs(ef - ec) / abs(ec))
        worst_m = max(worst_m, abs(nf - nc) / abs(nc))
    log(f'{name:18s}: max relative difference of the energy a(.,.) after prolongation = {worst_a:.3e}, '
        f'of the L2 norm^2 = {worst_m:.3e}')
out.close()
