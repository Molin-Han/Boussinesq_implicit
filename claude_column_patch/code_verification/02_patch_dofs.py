"""
Items 2 and 3: what is inside a star patch and a vanka patch, in serial and in parallel?

Two modes
---------
--mode solver  (default)
    The REAL solver script is executed (exec of its source in this __main__, so that
    '__main__.HDivSchurPC' resolves) for one time step on a small mesh, with
    ASMStarPC.get_patches / ASMVankaPC.get_patches wrapped so that every call is recorded:
    the function space handed to the smoother by PCMG/PCFIELDSPLIT and the index sets
    returned by Firedrake's own patch code.  Nothing of the patch construction is
    re-implemented here.
--mode direct
    The mesh hierarchy and the spaces are built exactly as in the solver script and
    ASMStarPC().get_patches(V) / ASMVankaPC().get_patches(V) are called directly on
    V = V_2D * Vy * Vb of every level with smoother (no solve, cheap, used for nz = 20).

What is reported, per MG level and per rank
-------------------------------------------
* number of index sets returned, how many are empty, sizes of the non-empty ones;
* every DoF of every patch is classified by field (u, u_y, b), by the entity it sits on
  (vertical facet / horizontal facet / cell interior) and located (x, z);
* every DoF gets a partition-independent key (field, kind of base entity, x of the base
  entity, position inside the column), so that patches of different runs can be compared;
* the patches are compared with the patches 'by definition'
      vanka(c) = DoFs of base cell c + DoFs of its two base vertices,
      star(v)  = DoFs of base vertex v + DoFs of the two base cells that contain v,
  built from the keys alone (no plex traversal), for every base cell / vertex of the
  global mesh: missing, duplicated and incomplete patches are counted;
* multiplicity of the covering (in how many patches does each DoF appear);
* which rank owns which base cell (the partition).

Usage (from the output folder; every run gets its own directory runs/patch_dofs_<tag>/):
    python 02_patch_dofs.py --patch vanka --nx 4 --nz 3 --support_check
    mpiexec --bind-to none -n 2 python 02_patch_dofs.py --patch vanka --nx 4 --nz 3
    python 02_patch_dofs.py --mode direct --patch star --nx 20 --nz 20
Results: runs/patch_dofs_<tag>/report.txt and patches.json (keys of all patches).
"""
import os
import sys
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--mode', choices=['solver', 'direct'], default='solver')
parser.add_argument('--patch', choices=['star', 'vanka'], default='vanka')
parser.add_argument('--nx', type=int, default=4)
parser.add_argument('--nz', type=int, default=3)
parser.add_argument('--degree', type=int, default=2)
parser.add_argument('--refinement', type=int, default=2)
parser.add_argument('--length', type=float, default=300000.0)
parser.add_argument('--height', type=float, default=4000.0)
parser.add_argument('--dt', type=float, default=50.0)
parser.add_argument('--shift', type=float, default=1.096622711232e-08)
parser.add_argument('--tag', type=str, default=None)
parser.add_argument('--support_check', action='store_true',
                    help='(serial) additionally check where the velocity functions of one patch are non-zero.')
my_args, _ = parser.parse_known_args()

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = "/home/mh2520/molin/phd/yr2/Boussinesq_implicit"
SOLVER = os.path.join(REPO, "Nonlinear_Boussinesq_Irk_SC.py")
sys.path.insert(0, REPO)
sys.path.insert(0, HERE)

import numpy as np                                           # noqa: E402
from mpi4py import MPI                                       # noqa: E402
import firedrake as fd                                       # noqa: E402
from firedrake.preconditioners import asm                    # noqa: E402
from firedrake.petsc import PETSc                            # noqa: E402
import utils                                                 # noqa: E402
import patch_tools as pt                                     # noqa: E402

_comm = MPI.COMM_WORLD
_rank, _size = _comm.rank, _comm.size
_L = float(my_args.length)
tag = my_args.tag or (f"{my_args.mode}_{my_args.patch}_deg{my_args.degree}_nx{my_args.nx}_nz{my_args.nz}_np{_size}")
RUN = os.path.join(HERE, "runs", "patch_dofs_" + tag)
if _rank == 0:
    os.makedirs(RUN, exist_ok=True)
_comm.Barrier()
os.chdir(RUN)

records = []       # one entry per call of get_patches


def record(pcname, V, ises):
    T = pt.describe_space(V, _L)
    mesh = T['mesh']
    dm = mesh.topology_dm
    patches = []
    for k, iset in enumerate(ises):
        idx = iset.indices.copy()
        if len(idx) == 0:
            patches.append(dict(empty=True, n=0))
            continue
        assert len(np.unique(idx)) == len(idx), "a DoF appears twice inside one patch"
        f, kd = T['field'][idx], T['kind'][idx]
        counts = {}
        for i in range(3):
            for kk in range(3):
                n = int(((f == i) & (kd == kk)).sum())
                if n:
                    counts[f"{i}:{kk}"] = n
        zs = T['z'][idx]
        bps = np.unique(T['base'][idx])
        patches.append(dict(
            empty=False, n=int(len(idx)), counts=counts,
            base_cells=sorted(pt.canon(T['px'][p], _L) for p in bps if T['iscell'][p]),
            base_vertices=sorted(pt.canon(T['px'][p], _L) for p in bps if not T['iscell'][p]),
            zmin=round(float(zs.min()), 6), zmax=round(float(zs.max()), 6),
            x_facet=sorted(set(pt.canon(v, _L) for v in T['x'][idx][(f == 0) & (kd == 0)])),
            x_nodes=sorted(set(pt.canon(v, _L) for v in T['x'][idx])),
            z_by_kind={f"{i}:{kk}": [float(T['z'][idx][(f == i) & (kd == kk)].min()),
                                     float(T['z'][idx][(f == i) & (kd == kk)].max()),
                                     int(len(np.unique(np.round(T['z'][idx][(f == i) & (kd == kk)], 6))))]
                       for i in range(3) for kk in range(3) if ((f == i) & (kd == kk)).any()},
            keys=sorted(pt.keys_of(T, idx)),
            idx=idx.tolist(),
        ))
    owned_idx = np.nonzero(T['owned'])[0]
    assert len(owned_idx) == sum(W.dof_dset.size * W.block_size for W in V)
    ps, pe = dm.getChart()
    cs, ce = dm.getHeightStratum(0)
    vs, ve = dm.getDepthStratum(0)
    nowned_cells = int((~T['ghost'][cs - ps:ce - ps]).sum())
    rec = dict(pc=pcname, rank=_rank,
               nlocal_cells=int(ce - cs), nlocal_vertices=int(ve - vs),
               nghost_cells=int(T['ghost'][cs - ps:ce - ps].sum()),
               nghost_vertices=int(T['ghost'][vs - ps:ve - ps].sum()),
               owned_cells=sorted(pt.canon(T['px'][c - ps], _L) for c in range(cs, ce) if not T['ghost'][c - ps]),
               owned_vertices=sorted(pt.canon(T['px'][v - ps], _L) for v in range(vs, ve) if not T['ghost'][v - ps]),
               nises=len(ises), patches=patches,
               owned_keys=pt.keys_of(T, owned_idx), owned_kind=T['kind'][owned_idx].tolist(),
               dx=float(_L / _comm.allreduce(nowned_cells, op=MPI.SUM)),
               T=T, V=V)
    records.append(rec)


# ----------------------------------------------------------------------------------------
# hooks around Firedrake's own patch construction
# ----------------------------------------------------------------------------------------
_orig_star = asm.ASMStarPC.get_patches
_orig_vanka = asm.ASMVankaPC.get_patches


def star_hook(self, V):
    ises = _orig_star(self, V)
    record('star', V, ises)
    return ises


def vanka_hook(self, V):
    ises = _orig_vanka(self, V)
    record('vanka', V, ises)
    return ises


asm.ASMStarPC.get_patches = star_hook
asm.ASMVankaPC.get_patches = vanka_hook

if my_args.mode == 'solver':
    sys.argv = [SOLVER,
                '--nx', str(my_args.nx), '--nz', str(my_args.nz), '--degree', str(my_args.degree),
                '--refinement', str(my_args.refinement),
                '--length', repr(my_args.length), '--height', repr(my_args.height),
                '--dt', repr(my_args.dt), '--tmax', repr(my_args.dt),          # one time step
                '--shift', repr(my_args.shift), '--maxit', '100', '--no_rotation', '--patch', my_args.patch,
                '--timing']                                                    # no VTK output, no monitors
    src = open(SOLVER).read()
    # The solver script does `from firedrake import *` (which exports e.g. ufl.rank) and defines many
    # globals of its own, so the names of this file are saved and restored around the exec.
    _saved = dict(globals())
    exec(compile(src, SOLVER, 'exec'), globals())
    globals().pop('print', None)                  # the solver script rebinds print to PETSc.Sys.Print
    globals().update(_saved)
else:
    nx, nz, deg = my_args.nx, my_args.nz, my_args.degree
    distribution_parameters = {"partition": True, "overlap_type": (fd.DistributedMeshOverlapType.VERTEX, 1)}
    m = fd.PeriodicIntervalMesh(nx, my_args.length, distribution_parameters=distribution_parameters)
    mh = fd.MeshHierarchy(m, refinement_levels=my_args.refinement)
    hierarchy = fd.ExtrudedMeshHierarchy(mh, my_args.height, layers=[nz] * (my_args.refinement + 1),
                                         extrusion_type='uniform')
    new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
    opts = PETSc.Options()
    for lev in range(1, my_args.refinement + 1):
        mesh = new_mh[lev]
        V = utils.extrude_RT(mesh, k=deg) * fd.FunctionSpace(mesh, 'DG', deg - 1) * utils.W_theta(mesh, k=deg)
        if my_args.patch == 'star':
            pc = asm.ASMStarPC()
            pc.prefix = f"direct_{lev}_pc_star_"
            opts[pc.prefix + "construct_dim"] = 0
        else:
            pc = asm.ASMVankaPC()
            pc.prefix = f"direct_{lev}_pc_vanka_"
            opts[pc.prefix + "construct_dim"] = 1
        pc.get_patches(V)

# ----------------------------------------------------------------------------------------
# analysis
# ----------------------------------------------------------------------------------------
lines = []


def out(s=""):
    lines.append(s)


def strip(rec):
    return {k: v for k, v in rec.items() if k not in ('T', 'V')}


all_recs = _comm.gather([strip(r) for r in records], root=0)
ncalls = len(records)

if _rank == 0:
    out(f"mode = {my_args.mode}, patch = {my_args.patch}, degree = {my_args.degree}, nx = {my_args.nx} "
        f"(x {2**my_args.refinement} on the finest level), nz = {my_args.nz}, MPI ranks = {_size}")
    out(f"number of calls of get_patches per rank: {[len(r) for r in all_recs]}")
    k = my_args.degree
    nz = my_args.nz
    expect = {'vanka': nz * (4 * k * k + k) + 2 * k, 'star': nz * (8 * k * k - k) + 4 * k}
    out(f"expected patch size for degree {k}, nz {nz}: vanka nz(4k^2+k)+2k = {expect['vanka']}, "
        f"star nz(8k^2-k)+4k = {expect['star']}" + ("   (= 18 nz + 4, 30 nz + 8)" if k == 2 else ""))
    summary = []
    for call in range(ncalls):
        recs = [all_recs[r][call] for r in range(_size)]
        pcname = recs[0]['pc']
        dx = recs[0]['dx']
        ncols = int(round(my_args.length / dx))
        out("")
        out("=" * 100)
        out(f"call {call}: {pcname} patches on the level with {ncols} base cells (dx = {dx:g} m)")
        tot_ises = sum(r['nises'] for r in recs)
        nonempty = [p for r in recs for p in r['patches'] if not p['empty']]
        nempty = tot_ises - len(nonempty)
        for r in recs:
            ne = sum(1 for p in r['patches'] if p['empty'])
            oc = r['owned_cells']
            out(f"  rank {r['rank']}: local base cells {r['nlocal_cells']} (ghost {r['nghost_cells']}), "
                f"local base vertices {r['nlocal_vertices']} (ghost {r['nghost_vertices']}), "
                f"index sets returned {r['nises']}, empty {ne}, non-empty {r['nises'] - ne}")
            # the owned cells as intervals of x (the partition)
            idxs = sorted(int(round((xc - 0.5 * dx) / dx)) for xc in oc)
            runs_ = []
            for i in idxs:
                if runs_ and i == runs_[-1][1] + 1:
                    runs_[-1][1] = i
                else:
                    runs_.append([i, i])
            out(f"          owns the base cells number {', '.join(f'{a}..{b}' for a, b in runs_)} "
                f"(counted from x = 0), {len(oc)} cells, and {len(r['owned_vertices'])} base vertices")
        sizes = sorted(set(p['n'] for p in nonempty))
        out(f"  total index sets (= 'total subdomain blocks' of PCASM) {tot_ises}, empty {nempty}, "
            f"non-empty {len(nonempty)}, sizes of the non-empty ones {sizes}")

        # patches by definition, from the keys of all owned DoFs
        allkeys = [tuple(kk) for r in recs for kk in r['owned_keys']]
        allkind = [kk for r in recs for kk in r['owned_kind']]
        assert len(set(allkeys)) == len(allkeys), "DoF keys are not unique"
        by_base = {}
        for kk in allkeys:
            by_base.setdefault((kk[1], kk[2]), []).append(kk)
        cells = sorted(b[1] for b in by_base if b[0] == 'c')
        verts = sorted(b[1] for b in by_base if b[0] == 'v')
        out(f"  global base cells {len(cells)}, global base vertices {len(verts)}, global DoFs {len(allkeys)}")
        expected = {}
        if pcname == 'vanka':
            for xc in cells:
                xl = pt.canon(xc - 0.5 * dx, _L)
                xr = pt.canon(xc + 0.5 * dx, _L)
                expected[('c', xc)] = frozenset(by_base[('c', xc)] + by_base[('v', xl)] + by_base[('v', xr)])
        else:
            for xv in verts:
                xl = pt.canon(xv - 0.5 * dx, _L)
                xr = pt.canon(xv + 0.5 * dx, _L)
                expected[('v', xv)] = frozenset(by_base[('v', xv)] + by_base[('c', xl)] + by_base[('c', xr)])
        got = [frozenset(tuple(kk) for kk in p['keys']) for p in nonempty]
        exp_sets = set(expected.values())
        n_match = sum(1 for g in got if g in exp_sets)
        n_dupl = len(got) - len(set(got))
        n_missing = len(exp_sets - set(got))
        n_wrong = sum(1 for g in got if g not in exp_sets)
        out(f"  comparison with the patches by definition: expected {len(exp_sets)}, found identical {n_match}, "
            f"missing {n_missing}, not among the expected (incomplete or wrong) {n_wrong}, duplicated {n_dupl}")
        mult = {}
        for g in got:
            for kk in g:
                mult[kk] = mult.get(kk, 0) + 1
        kind_of = dict(zip(allkeys, allkind))
        table = {}
        for kk in allkeys:
            name = (pt.FIELD[kk[0]], pt.KIND[kind_of[kk]])
            table.setdefault(name, set()).add(mult.get(kk, 0))
        out("  covering multiplicity (number of patches that contain a DoF):")
        for name in sorted(table):
            out(f"      {name[0]:18s} {name[1]:18s}: {sorted(table[name])}")
        summary.append(dict(call=call, pc=pcname, ncols=ncols, tot_ises=tot_ises, empty=nempty,
                            nonempty=len(nonempty), sizes=sizes, expected=len(exp_sets), identical=n_match,
                            missing=n_missing, wrong=n_wrong, duplicated=n_dupl,
                            owned_cells_per_rank=[len(r['owned_cells']) for r in recs]))

        # one patch in detail (the first non-empty patch of rank 0 that does not touch the periodic seam)
        def span(p):
            xs_ = p['base_cells'] + p['base_vertices']
            return max(xs_) - min(xs_)
        cand = [p for p in recs[0]['patches'] if not p['empty'] and span(p) < 3 * dx]
        p = cand[0] if cand else nonempty[0]
        out(f"  one {pcname} patch in detail: {p['n']} DoFs")
        out(f"      base cells (x of the midpoint)  : {p['base_cells']}")
        out(f"      base vertices (x)               : {p['base_vertices']}")
        out(f"      x of the vertical-facet DoFs    : {p['x_facet']}")
        out(f"      x of all nodes                  : {p['x_nodes']}")
        out(f"      z range of all DoFs             : [{p['zmin']:g}, {p['zmax']:g}]  (L_z = {my_args.height:g})")
        out(f"      {'field':18s} {'entity':18s} {'DoFs':>6s}   z_min      z_max   distinct z")
        for key in sorted(p['counts']):
            i, kk = (int(s) for s in key.split(':'))
            zmin, zmax, nzd = p['z_by_kind'][key]
            out(f"      {pt.FIELD[i]:18s} {pt.KIND[kk]:18s} {p['counts'][key]:6d}   {zmin:8.2f}  {zmax:9.2f}   {nzd}")
    out("")
    out("SUMMARY " + json.dumps(summary))
    with open("patches.json", "w") as f:
        json.dump([[[{kk: v for kk, v in p.items() if kk != 'idx'} for p in all_recs[r][c]['patches']]
                    for r in range(_size)] for c in range(ncalls)], f)

# ----------------------------------------------------------------------------------------
# optional: where are the velocity functions of one patch non-zero?  (serial only)
# ----------------------------------------------------------------------------------------
if my_args.support_check and _size == 1 and ncalls > 0:
    rec = records[-1]                      # finest level
    T, V = rec['T'], rec['V']
    mesh = T['mesh']
    dx = rec['dx']
    cand = [p for p in rec['patches'] if not p['empty']
            and max(p['base_cells'] + p['base_vertices']) - min(p['base_cells'] + p['base_vertices']) < 3 * dx]
    p = cand[0]
    idx = np.array(p['idx'])
    V0 = V[0]
    DG0 = fd.FunctionSpace(mesh, "DG", 0)
    xcell = fd.Function(DG0).interpolate(fd.SpatialCoordinate(mesh)[0]).dat.data_ro
    q = fd.TestFunction(DG0)
    loc0 = V.dof_dset.local_ises[0].indices
    inv0 = {int(g): n for n, g in enumerate(loc0)}
    out("")
    out("=" * 100)
    out(f"support of the velocity functions of the detailed {rec['pc']} patch (finest level), serial")
    out(f"   patch: base cells {p['base_cells']}, base vertices {p['base_vertices']}")
    for kk in range(3):
        sel = idx[(T['field'][idx] == 0) & (T['kind'][idx] == kk)]
        if len(sel) == 0:
            continue
        u = fd.Function(V0)
        rng = np.random.default_rng(1)
        u.dat.data[[inv0[int(g)] for g in sel]] = 1.0 + rng.random(len(sel))
        ex = fd.assemble(q * u[0] ** 2 * fd.dx).dat.data_ro
        ey = fd.assemble(q * u[1] ** 2 * fd.dx).dat.data_ro
        ez = fd.assemble(q * u[2] ** 2 * fd.dx).dat.data_ro
        dv = fd.assemble(q * fd.div(u) ** 2 * fd.dx).dat.data_ro
        scale = max(ex.max(), ez.max())             # int_cell u_x^2, int_cell u_z^2 on the same scale
        cols_x = sorted(set(pt.canon(v, _L) for v in xcell[ex > 1e-12 * scale]))
        cols_z = sorted(set(pt.canon(v, _L) for v in xcell[ez > 1e-12 * scale]))
        cols_d = sorted(set(pt.canon(v, _L) for v in xcell[dv > 1e-12 * dv.max()])) if dv.max() > 0 else []
        out(f"   velocity DoFs on {pt.KIND[kk]:16s} ({len(sel):4d} DoFs): u_x non-zero in the columns {cols_x}, "
            f"u_z non-zero in the columns {cols_z}, max |u_y| {np.sqrt(ey.max()):.1e}; "
            f"div u non-zero in the columns {cols_d}")

if _rank == 0:
    text = "\n".join(lines)
    sys.stdout.write(text + "\n")
    sys.stdout.flush()
    with open("report.txt", "w") as f:
        f.write(text + "\n")
