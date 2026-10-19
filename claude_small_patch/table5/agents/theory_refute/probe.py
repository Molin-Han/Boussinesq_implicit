# theory_refute probe: FS vs MONO shifted preconditioners of LB_time_slice.py (Table 5 settings),
# small serial meshes only.  The problem set-up below is copied verbatim (same forms) from
# /home/mh2520/molin/phd/yr2/Boussinesq_implicit/LB_time_slice.py (lines 98-212, 214-398).
import sys, os, time, argparse, re
import numpy as np
REPO = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit'
sys.path.insert(0, REPO)
from firedrake import *
from firedrake.petsc import PETSc
from firedrake import dmhooks
from firedrake.utils import unique
from contextlib import ExitStack
from itertools import chain
import utils

print = PETSc.Sys.Print
ap = argparse.ArgumentParser()
ap.add_argument('--nx', type=int, default=8)
ap.add_argument('--nz', type=int, default=8)
ap.add_argument('--solvers', type=str, default='fs,mono')
ap.add_argument('--eigs', type=str, default=None, help='fixed Chebyshev "emin,emax" on all MG levels (both solvers)')
ap.add_argument('--pc_side', type=str, default='left')
ap.add_argument('--steps', type=int, default=2)
ap.add_argument('--probe', action='store_true')
ap.add_argument('--tag', type=str, default='run')
ap.add_argument('--mono_use_amat', type=str, default=None, help='0/1: set pc_use_amat on the MONO inner PCMG (PETSc default for PCMG is 1)')
ap.add_argument('--fs_use_amat', type=str, default=None, help='0/1: set pc_use_amat on the FS inner PCMG')
args, _ = ap.parse_known_args()

nx, nz = args.nx, args.nz
length, height, deg = 300000.0, 4000.0, 2
dt = Constant(100.0); shift = Constant(1e-4); delta = 1e-4
use_rotation = False
appctx = {"dt": dt, "shift": shift}
U_MEAN = 0.0

def vector_3D(u, uy):
    return u + uy * utils.j()

class HDivSchurPC(AuxiliaryOperatorPC):
    _prefix = "helmholtzschurpc_"
    def form(self, pc, v, u):
        prefix = (pc.getOptionsPrefix() or "") + self._prefix
        rotation = PETSc.Options().getBool(f"{prefix}use_rotation", False)
        appctx_PC = self.get_appctx(pc)
        dtc = appctx_PC["dt"]; dlt = appctx_PC["shift"]; n = appctx_PC["n"]
        One = as_vector([1., 1., 1.])
        uxz, uy, b = split(u)
        wxz, wy, q = split(v)
        velo = vector_3D(uxz, uy)
        unph = Constant(0.5) * (velo + One)
        bnph = Constant(0.5) * (b + Constant(1.))
        w = vector_3D(wxz, wy)
        pnp1 = - Constant(1.) / dlt * div(velo)
        Jp_ = lhs(utils.LB_velocity(velo, One, unph, w, bnph, pnp1, n, dtc, use_rotation=rotation, twoD=False, U_mean=Constant(U_MEAN)))
        Jp_ += lhs(utils.LB_buoyancy(b, Constant(1.), q, unph, bnph, n, dtc, twoD=False, U_mean=Constant(U_MEAN)))
        _, bcs_ = super().form(pc, u, v)
        return (Jp_, bcs_)

class MonolithicShiftPC(AuxiliaryOperatorPC):
    _prefix = "monoshiftpc_"
    def form(self, pc, v, u):
        test, trial = Jp.arguments()
        a = replace(Jp, {test: v, trial: u})
        _, bcs_ = super().form(pc, u, v)
        return (a, bcs_)

class MonoBackSubPC(MonolithicShiftPC):
    """MONO one cycle, then the pressure is REPLACED by the exact back-substitution
    x_p = (delta M)^{-1}(r_p - D x_v)  (what PCFIELDSPLIT 'full' does)."""
    def initialize(self, pc):
        super().initialize(pc)
        Pm = self.P.petscmat
        self.isv = PETSc.IS().createGeneral(np.concatenate([W._ises[i].indices for i in range(3)]).astype(PETSc.IntType), comm=PETSc.COMM_SELF)
        self.isp = W._ises[3]
        self.D = Pm.createSubMatrix(self.isp, self.isv)
        self.Mpp = Pm.createSubMatrix(self.isp, self.isp)
        self.kspM = PETSc.KSP().create(comm=PETSc.COMM_SELF)
        self.kspM.setOperators(self.Mpp); self.kspM.setType('preonly'); self.kspM.getPC().setType('lu'); self.kspM.setUp()
    def apply(self, pc, x, y):
        super().apply(pc, x, y)
        xp = x.getSubVector(self.isp); yv = y.getSubVector(self.isv)
        t = xp.duplicate(); self.D.mult(yv, t); t.aypx(-1.0, xp)
        x.restoreSubVector(self.isp, xp); y.restoreSubVector(self.isv, yv)
        yp = y.getSubVector(self.isp); self.kspM.solve(t, yp); y.restoreSubVector(self.isp, yp)

distribution_parameters = {"partition": True, "overlap_type": (DistributedMeshOverlapType.VERTEX, 1)}
m = PeriodicIntervalMesh(nx, length, distribution_parameters=distribution_parameters)
mh = MeshHierarchy(m, refinement_levels=2)
hierarchy = ExtrudedMeshHierarchy(mh, height, layers=[nz] * 3, extrusion_type='uniform')
new_mh = utils.high_dim_mesh_hierarchy(hierarchy, dim=3)
mesh = new_mh[-1]
mesh.name = "finest"
x, y, z = SpatialCoordinate(mesh)
V_2D = utils.extrude_RT(mesh, k=deg)
Vy = FunctionSpace(mesh, 'DG', deg-1)
Pressure = FunctionSpace(mesh, 'DG', deg-1)
Vb = utils.W_theta(mesh, k=deg)
W = V_2D * Vy * Vb * Pressure
Un = Function(W); Unp1 = Function(W)
unxz, uny, bn, pn = split(Un)
unp1xz, unp1y, bnp1, pnp1 = split(Unp1)
w_xz, wy, q, phi = TestFunctions(W)
xc = Constant(length/2); a = Constant(5000)
u0_slice, u0yic, b0ic, p0ic = Un.subfunctions
b0ic.project(3e-4*sin(pi*z/height)/(1+((x-xc)**2)/a**2))
IC = Function(W); IC.assign(Un)
bc1 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "top")
bc2 = DirichletBC(W.sub(0), as_vector([0., 0., 0.]), "bottom")
bcs = [bc1, bc2]
un = vector_3D(unxz, uny); unp1 = vector_3D(unp1xz, unp1y)
unph = 0.5 * (un+unp1); bnph = 0.5 * (bn+bnp1)
w = vector_3D(w_xz, wy)
n = FacetNormal(mesh)
appctx.update({"n": n})
eqn = utils.LB_velocity(unp1, un, unph, w, bnph, pnp1, n, dt, use_rotation=use_rotation, U_mean=Constant(U_MEAN))
eqn += utils.LB_buoyancy(bnp1, bn, q, unph, bnph, n, dt, U_mean=Constant(U_MEAN))
eqn += utils.LB_pressure(unp1, phi)
shift_eqn = eqn + shift * pnp1 * phi * dx
Jp = derivative(shift_eqn, Unp1)
v_basis = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), W.sub(1), W.sub(2), v_basis])

patch_params = {"pc_type": "python", "pc_python_type": "firedrake.ASMStarPC", "pc_star_construct_dim": 0,
                "pc_star_sub_sub_pc_type": "lu",
                "pc_star_sub_sub_pc_factor_mat_ordering_type": "rcm", "pc_star_sub_sub_pc_factor_reuse_ordering": None}
cheb = {'ksp_type': 'chebyshev', 'ksp_max_it': 6}
if args.eigs is not None:
    cheb['ksp_chebyshev_eigenvalues'] = args.eigs
helm = {'pc_type': 'mg', 'pc_mg_type': 'full', 'pc_mg_cycle_type': 'v',
        'mg_levels': {**patch_params, **cheb},
        'mg_coarse': {'ksp_type': 'preonly', 'pc_type': 'lu'}}
base = {'snes_type': 'ksponly', 'ksp_type': 'gmres', 'ksp_atol': 0, 'ksp_rtol': 1e-8, 'ksp_max_it': 100,
        'ksp_converged_maxits': None, 'ksp_pc_side': args.pc_side}
if args.fs_use_amat is not None:
    helm['pc_use_amat'] = bool(int(args.fs_use_amat))
P_FS = {**base, 'pc_type': 'fieldsplit', 'pc_fieldsplit_type': 'schur', 'pc_fieldsplit_schur_fact_type': 'full',
        'pc_fieldsplit_0_fields': '3', 'pc_fieldsplit_1_fields': '0,1,2',
        'fieldsplit_0': {'ksp_type': 'preonly', 'pc_type': 'bjacobi', 'sub_pc_type': 'ilu'},
        'fieldsplit_1': {'use_rotation': use_rotation, 'ksp_type': 'preonly', 'pc_type': 'python',
                         'pc_python_type': __name__ + '.HDivSchurPC', 'helmholtzschurpc': helm}}
mono_inner = {'mat_type': 'aij', 'pc_type': 'mg', 'pc_mg_type': 'full', 'pc_mg_cycle_type': 'v',
              'mg_levels': {**patch_params, **cheb}, 'mg_coarse': {'ksp_type': 'preonly', 'pc_type': 'lu'}}
if args.mono_use_amat is not None:
    mono_inner['pc_use_amat'] = bool(int(args.mono_use_amat))
P_MONO = {**base, 'pc_type': 'python', 'pc_python_type': __name__ + '.MonolithicShiftPC', 'monoshiftpc': mono_inner}
P_MONOBS = {**base, 'pc_type': 'python', 'pc_python_type': __name__ + '.MonoBackSubPC', 'monoshiftpc': mono_inner}
CONF = {'fs': P_FS, 'mono': P_MONO, 'monobs': P_MONOBS}

nprob = NonlinearVariationalProblem(eqn, Unp1, bcs=bcs, Jp=Jp)

def solver_ctx(s):
    problem = s._problem
    forms = (problem.F, problem.J, problem.Jp)
    coefficients = unique(chain.from_iterable(f.coefficients() for f in forms if f is not None))
    sdm = s.snes.getDM(); dms = []
    for c in chain(coefficients, problem.dirichlet_bcs()):
        dm_ = c.function_space().dm
        if dm_ == sdm: continue
        if dm_ not in dms: dms.append(dm_)
    dms.append(sdm)
    st = ExitStack()
    for c in chain([s.inserted_options()], [dmhooks.add_hooks(d, s, appctx=s._ctx) for d in dms], s._transfer_operators):
        st.enter_context(c)
    return st

def inner_mg(s, name):
    pc = s.snes.ksp.getPC()
    if name == 'fs':
        k1 = pc.getFieldSplitSubKSP()[1]
        return k1.getPC().getPythonContext().pc
    return pc.getPythonContext().pc

def cheb_info(s, name):
    mg = inner_mg(s, name)
    out = []
    for l in range(1, mg.getMGLevels()):
        ksp = mg.getMGSmoother(l)
        fn = os.path.join(os.getcwd(), f'_view_{args.tag}_{name}_{l}.txt')
        vw = PETSc.Viewer().createASCII(fn, 'w', comm=PETSc.COMM_SELF); ksp.view(vw); vw.destroy()
        txt = open(fn).read(); os.remove(fn)
        tg = re.search(r'eigenvalue targets used: min ([-0-9.e+]+), max ([-0-9.e+]+)', txt)
        es = re.search(r'eigenvalues estimated via \w+: min ([-0-9.e+]+), max ([-0-9.e+]+)', txt)
        out.append((l, tg.groups() if tg else None, es.groups() if es else None))
    return out

solvers = {}
res = {}
for name in args.solvers.split(','):
    t0 = time.time()
    s = NonlinearVariationalSolver(nprob, nullspace=nullspace, solver_parameters=CONF[name], appctx=appctx,
                                   options_prefix=f'{name}_')
    its = []
    Un.assign(IC); Unp1.assign(Un)
    for step in range(args.steps):
        s.solve()
        its.append(s.snes.ksp.getIterationNumber())
        if step == 0:
            sol1 = Function(W); sol1.assign(Unp1)
        Un.assign(Unp1)
    solvers[name] = (s, sol1)
    res[name] = its
    print(f'[{args.tag}] nx={nx} nz={nz} solver={name} pc_side={args.pc_side} eigs={args.eigs} mono_use_amat={args.mono_use_amat} fs_use_amat={args.fs_use_amat} its/step={its} wall={time.time()-t0:.1f}s')
    mgx = inner_mg(s, 'fs' if name == 'fs' else 'mono'); kf = mgx.getMGSmoother(mgx.getMGLevels()-1); Afine, Pfine = kf.getOperators()
    print(f'[{args.tag}]   {name} finest-level smoother: Amat is outer J: {Afine.handle == s.snes.ksp.getOperators()[0].handle}; Amat type {Afine.getType()}; Amat==Pmat: {Afine.handle == Pfine.handle}')
    if name in ('fs', 'mono', 'monobs'):
        for l, tg, es in cheb_info(s, 'fs' if name == 'fs' else 'mono'):
            print(f'[{args.tag}]   {name} MG level {l}: cheb targets {tg}  estimated(min,max) {es}')

if args.probe and 'fs' in solvers and 'mono' in solvers:
    sF = solvers['fs'][0]; sM = solvers['mono'][0]
    Jmat = sF.snes.ksp.getOperators()[0]
    ivs = [W._ises[i].indices for i in range(4)]
    iv = np.concatenate(ivs[:3]); ip = ivs[3]
    bcn = np.concatenate([W._ises[0].indices[bc.nodes] for bc in bcs])
    rng = np.random.RandomState(1234)
    # pressure mass for (delta M)^{-1}
    pp, qq = TrialFunction(Pressure), TestFunction(Pressure)
    Mp = assemble(pp*qq*dx).petscmat
    kM = PETSc.KSP().create(comm=PETSc.COMM_SELF); kM.setOperators(Mp); kM.setType('preonly'); kM.getPC().setType('lu'); kM.setUp()
    def dMinv(vec_np):
        b_ = Mp.createVecLeft(); b_.setArray(vec_np); x_ = Mp.createVecRight(); kM.solve(b_, x_); return x_.getArray().copy() / delta
    def applyP(s, rnp):
        xv = Jmat.createVecLeft(); xv.setArray(rnp); yv = Jmat.createVecRight()
        with solver_ctx(s):
            s.snes.ksp.getPC().apply(xv, yv)
        out = yv.getArray().copy()
        out[ip] -= out[ip].mean()
        return out
    def mk(kind):
        r = rng.randn(Jmat.getSize()[0])
        r[bcn] = 0.0
        if kind == 'v': r[ip] = 0.0
        if kind == 'p': r[iv] = 0.0
        return r
    nrm = np.linalg.norm
    zs = {}
    for kind in ('v', 'p', 'vp'):
        r = mk(kind)
        yF = applyP(sF, r); yM = applyP(sM, r)
        dv = nrm(yM[iv]-yF[iv])/nrm(yF[iv]); dp = nrm(yM[ip]-yF[ip])/nrm(yF[ip])
        print(f'[{args.tag}] probe r={kind}: |xv_M-xv_F|/|xv_F|={dv:.3e}  |xp_M-xp_F|/|xp_F|={dp:.3e}')
        zs[kind] = (r, yM[ip]-yF[ip])
    # z independence of r_v: z(r_v + r_p) vs z(r_p) for the SAME r_p
    r_p_only = zs['p'][0].copy(); r_v_only = zs['v'][0].copy()
    r_sum = r_p_only + r_v_only
    zsum = applyP(sM, r_sum)[ip] - applyP(sF, r_sum)[ip]
    zp = zs['p'][1]
    print(f'[{args.tag}] probe linearity/independence: |z(rv+rp) - z(rp)|/|z(rp)| = {nrm(zsum - zp)/nrm(zp):.3e}; |z(rv)|/|xp_F(rv)| = {nrm(zs["v"][1])/nrm(applyP(sF, r_v_only)[ip]):.3e}')
    dm_rp = dMinv(r_p_only[ip]); dm_rp -= dm_rp.mean()
    print(f'[{args.tag}] probe |z(rp)|/|(dM)^-1 rp| = {nrm(zp)/nrm(dm_rp):.3e}; |xp_F(rp)|/|(dM)^-1 rp| = {nrm(applyP(sF, r_p_only)[ip])/nrm(dm_rp):.3e}')
    # z on the actual Krylov direction: r = A_0 * (P_F^{-1} b) with b the first-step RHS
    np.save(os.path.join(os.getcwd(), f'z_{args.tag}.npy'), zp)
