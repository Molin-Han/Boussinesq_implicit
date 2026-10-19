'''
Verification of the exact modal reduction (modal1d.py) against the actual Firedrake
discretisation and the actual patches of firedrake.ASMStarPC (construct_dim 0) and
firedrake.ASMVankaPC (construct_dim 1).

1. Velocity-only operator a(u, w) = (u_x, w_x) + m (u_z, w_z) + c (div u, div w) on the
   space V_2D of the solver (utils.extrude_RT on the embedded extruded periodic mesh):
   the eigenvalues of B_ASM A computed from the assembled 2D matrix and the patches built
   by firedrake/preconditioners/asm.py are compared with the union over the vertical
   modes n of the eigenvalues of the 1D modal problems.
2. The full field-1 operator on (u, u_y, b) (stage-derivative form, linearised about rest,
   f = 0): its ASM-preconditioned spectrum is computed with the same patches, to show that
   the additional fields (u_y, b) do not change the small eigenvalues.

Run (serial):  python verify_modal_reduction_firedrake.py
Writes verify_modal_reduction_firedrake.out
'''
import sys
import os
import numpy as np
import scipy.linalg as sla

sys.path.insert(0, '/home/mh2520/molin/phd/yr2/Boussinesq_implicit')
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)

from firedrake import *                       # noqa
from firedrake.preconditioners.asm import (build_star_indices, build_vanka_indices,
                                           get_local_ises_indices)
import utils                                   # the solver's utils.py (read only)
import modal1d

out = open(os.path.join(here, 'verify_modal_reduction_firedrake.out'), 'w')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    out.write(s + '\n')
    out.flush()


# ---------------------------------------------------------------- parameters
Nc, Nz, deg = 8, 6, 2
dx_, Lz = 3750.0, 4000.0
dt, Nfreq, K = 50.0, 1.0e-2, 1.0
Lx_full = 3.0e5                                  # only used for C_M (the shift)
g = modal1d.groups(dt, Lx_full, Lz, int(round(Lx_full / dx_)), Nz, K)
c, m, delta = g['c'], g['m'], g['delta']
log(f'Nc = {Nc}, Nz = {Nz}, degree = {deg}, dx = {dx_}, Lz = {Lz}, c = {c:.6e}, m = {m}, '
    f'delta = {delta:.6e}, alpha = c/dx^2 = {c/dx_**2:.4f}')

# ---------------------------------------------------------------- mesh as in the solver
distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
base = PeriodicIntervalMesh(Nc, Nc * dx_, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(base, 0)
hierarchy = ExtrudedMeshHierarchy(mh, Lz, layers=[Nz], extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]
n = FacetNormal(mesh)

V = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg - 1)
Vb = utils.W_theta(mesh, k=deg)


def dense(a, bcs):
    A = assemble(a, bcs=bcs, mat_type='aij')
    P = A.petscmat
    ia, ja, va = P.getValuesCSR()
    import scipy.sparse as sp
    return sp.csr_matrix((va, ja, ia), shape=P.getSize()).toarray()


def patches_of(W, kind):
    dm = W.mesh().unique().topology_dm
    ises = get_local_ises_indices(W)
    out_ = []
    if kind == 'star':
        for seed in range(*dm.getDepthStratum(0)):
            iset = build_star_indices(W, ises, dm, 'natural', 'dummy_', seed)
            out_.append(np.array(iset.indices))
    else:
        Z = (tuple(W[i] for i in range(len(W))), tuple())
        Zi = (tuple(ises[i] for i in range(len(W))), tuple())
        for seed in range(*dm.getDepthStratum(1)):
            iset = build_vanka_indices(Z, Zi, dm, 'natural', 'dummy_', True, seed)
            out_.append(np.array(iset.indices))
    return out_


def bc_dofs(W, bcs):
    off = np.cumsum([0] + [Wi.dof_dset.layout_vec.getSize() for Wi in W])
    d = []
    for bc in bcs:
        idx = bc.function_space().index
        idx = 0 if idx is None else idx
        d += list(off[idx] + np.array(bc.nodes) * W[idx].block_size if W[idx].block_size == 1
                  else off[idx] + np.array(bc.nodes))
    return np.unique(np.array(d, dtype=int))


def asm_spectrum(A, patches, drop):
    B = np.zeros_like(A)
    for p in patches:
        B[np.ix_(p, p)] += np.linalg.inv(A[np.ix_(p, p)])
    keep = np.setdiff1d(np.arange(A.shape[0]), drop)
    BA = (B @ A)[np.ix_(keep, keep)]
    return np.sort_complex(np.linalg.eigvals(BA))


# ================================================================ 1. velocity only
u, w = TrialFunction(V), TestFunction(V)
kz = utils.k()
a_vel = (inner(u, w) + Constant(m - 1.0) * inner(u, kz) * inner(w, kz)
         + Constant(c) * div(u) * div(w)) * dx
bcs_v = [DirichletBC(V, as_vector([0., 0., 0.]), "top"),
         DirichletBC(V, as_vector([0., 0., 0.]), "bottom")]
A2 = dense(a_vel, bcs_v)
bcd = np.unique(np.concatenate([bc.nodes for bc in bcs_v]))
log(f'velocity space: {A2.shape[0]} dofs, {len(bcd)} boundary dofs, symmetric defect {np.abs(A2-A2.T).max():.2e}')

lam2 = modal1d.vertical_eigenvalues(Nz, Lz, deg)
L1 = modal1d.Level(Nc, dx_, deg)
res = {}
for kind in ('star', 'vanka'):
    pts = patches_of(V, kind)
    sizes = sorted(set(len(p) for p in pts))
    ev2 = asm_spectrum(A2, pts, bcd)
    # modal prediction
    ev1 = []
    A0 = L1.matrix(c, 0.0, m, with_w=False)
    B0 = modal1d.asm_matrix(A0, L1.patches(kind, with_w=False))
    ev1 += list(np.linalg.eigvals(B0 @ A0).real)
    for l2 in lam2:
        An = L1.matrix(c, np.sqrt(l2), m)
        Bn = modal1d.asm_matrix(An, L1.patches(kind))
        ev1 += list(np.linalg.eigvals(Bn @ An).real)
    ev1 = np.sort(np.array(ev1))
    ev2r = np.sort(ev2.real)
    err = np.abs(ev1 - ev2r).max()
    log(f'[velocity only] {kind:5s}: {len(pts)} patches of sizes {sizes}; '
        f'2D spectrum in [{ev2r[0]:.6e}, {ev2r[-1]:.6f}] ({len(ev2r)} eigenvalues, max |Im| = {np.abs(ev2.imag).max():.1e}); '
        f'modal 1D spectrum in [{ev1[0]:.6e}, {ev1[-1]:.6f}] ({len(ev1)}); '
        f'max |difference of sorted spectra| = {err:.3e}')
    res[kind] = ev2r

# ================================================================ 2. full field 1: (u, u_y, b)
W1 = V * Vy * Vb
(u, uy, b) = TrialFunctions(W1)
(w, wy, q) = TestFunctions(W1)
N2 = Constant(Nfreq ** 2)
hdt = Constant(dt / 2.0)
a_f1 = (inner(u, w) + uy * wy + hdt / Constant(delta) * div(u) * div(w)
        - hdt * inner(w, kz) * b
        + q * b + hdt * N2 * q * inner(u, kz)) * dx
bcs_1 = [DirichletBC(W1.sub(0), as_vector([0., 0., 0.]), "top"),
         DirichletBC(W1.sub(0), as_vector([0., 0., 0.]), "bottom")]
A3 = dense(a_f1, bcs_1)
log(f'field-1 space (u, u_y, b): {A3.shape[0]} dofs = {A3.shape[0]//Nc} per column')
for kind in ('star', 'vanka'):
    pts = patches_of(W1, kind)
    sizes = sorted(set(len(p) for p in pts))
    ev3 = asm_spectrum(A3, pts, bcd)          # velocity is field 0, so the bc dofs keep their numbers
    mod = np.abs(ev3)
    small = np.sort(ev3.real[mod < 0.5])
    small_v = res[kind][res[kind] < 0.5]
    log(f'[field 1      ] {kind:5s}: {len(pts)} patches of sizes {sizes} (solver at Nz = 20: star 608 = 30 Nz + 8, vanka 364 = 18 Nz + 4; '
        f'here 30 Nz + 8 = {30*Nz+8}, 18 Nz + 4 = {18*Nz+4}); '
        f'spectrum: min Re = {ev3.real.min():.6e}, max Re = {ev3.real.max():.6f}, max |Im| = {np.abs(ev3.imag).max():.3e}; '
        f'eigenvalues with |.| < 0.5: {len(small)} (velocity only: {len(small_v)})')
    if len(small) and len(small) == len(small_v):
        log(f'                 max relative difference of the small eigenvalues to the velocity-only ones: '
            f'{np.abs(small - small_v).max() / small_v.max():.3e}')
out.close()
