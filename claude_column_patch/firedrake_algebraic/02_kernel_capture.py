'''
02  Kernel capture of the star and of the column (vanka) patches.

For the in-plane velocity space V (RT_k on the extruded periodic mesh, u.n = 0 on top and
bottom) and the pressure space Q = DG_{k-1}, with B[i, j] = (phi_i, div u_j):

  * dim ker(B) on the free velocity DoFs (SVD of the dense matrix)
  * for every patch of the mixed space (u, u_y, b) built by the real Firedrake classes,
    the velocity DoFs of the patch and the kernel of B restricted to the fields spanned
    by them (all other DoFs zero)
  * the dimension of the span of all patch kernels (rank of the matrix whose columns are
    the patch-kernel basis vectors extended by zero)
  * the same for the "degraded" divergence d_z u_z (the premise of the hypothesis)
  * the dependence of B on the geometry (L_x, L_z)

Expected from the stream-function picture (derivation in README.md):
    dim ker(div)           = k ncol (k nz - 1) + 1
    star patch kernel      = (2k - 1)(k nz - 1)      span of all = k ncol (k nz - 1)
    column patch kernel    = (k - 1)(k nz - 1)       span of all = (k - 1) ncol (k nz - 1)

usage: python 02_kernel_capture.py
output: results_02_kernel_capture.csv, results_02_kernel_capture.txt
'''
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fa_common as fa
from fa_common import *          # noqa

out = []


def report(line):
    PETSc.Sys.Print(line)
    out.append(line)


def kernel_data(nx, nz, length, height, deg, refinement=2):
    mh = fa.build_hierarchy(nx, nz, length, height, refinement=refinement)
    mesh = mh[-1]
    V, Vy, Vb, Q, W, W1 = fa.build_spaces(mesh, deg)
    ncol = mesh.num_cells()
    B = fa.assemble_scipy(fa.divergence_form(V, Q)).toarray()
    # degraded divergence: only the vertical derivative of the vertical velocity
    u = TrialFunction(V)
    phi = TestFunction(Q)
    Bz = fa.assemble_scipy(phi * grad(u)[2, 2] * dx).toarray()
    Bx = fa.assemble_scipy(phi * grad(u)[0, 0] * dx).toarray()
    bcn = fa.bc_dofs_velocity(V)
    free = np.setdiff1d(np.arange(V.dim()), bcn)
    f1 = fa.field_offsets(W1)
    res = {"ncol": ncol, "nz": nz, "deg": deg, "length": length, "height": height,
           "dimV": V.dim(), "n_bc": len(bcn), "n_free": len(free), "dimQ": Q.dim(),
           "B": B, "free": free, "max|B - Bx - Bz|": abs(B - Bx - Bz).max()}
    Bf = B[:, free]
    Bzf = Bz[:, free]
    kdim, s = fa.null_space_dim(Bf)
    res["dim_ker"] = kdim
    res["rank_B"] = len(free) - kdim
    res["sv_gap"] = (s[res["rank_B"] - 1] / s[0], (s[res["rank_B"]] / s[0]) if res["rank_B"] < len(s) else 0.0)
    res["dim_ker_dz"] = fa.null_space_dim(Bzf)[0]
    res["formula_ker"] = deg * ncol * (deg * nz - 1) + 1
    Nglob = fa.null_space_basis(Bf)
    for kind in ["star", "vanka"]:
        patches = fa.patch_index_sets(W1, kind)
        sizes = sorted(set(len(p) for p in patches))
        # velocity part of the patches, in the numbering of the free velocity DoFs
        # (field 0 of W1 occupies the first dim(V) indices of the monolithic numbering)
        assert np.array_equal(f1[0], np.arange(V.dim()))
        pu = fa.restrict_patches(patches, free)
        pdims = []
        pdims_dz = []
        cols = []
        cols_dz = []
        for p in pu:
            Nj = fa.null_space_basis(Bf[:, p])
            pdims.append(Nj.shape[1])
            Z = np.zeros((len(free), Nj.shape[1]))
            Z[p, :] = Nj
            cols.append(Z)
            Njz = fa.null_space_basis(Bzf[:, p])
            pdims_dz.append(Njz.shape[1])
            Zz = np.zeros((len(free), Njz.shape[1]))
            Zz[p, :] = Njz
            cols_dz.append(Zz)
        Z = np.hstack(cols) if cols else np.zeros((len(free), 0))
        if Z.shape[1] > 0:
            sz = np.linalg.svd(Z, compute_uv=False)
            span = int(np.sum(sz > 1e-10 * sz[0]))
            # how much of the global kernel is NOT in the span: project the orthonormal
            # global kernel basis on the orthogonal complement of the span
            Qz, _ = np.linalg.qr(Z[:, :])
            Uz, sz2, _ = np.linalg.svd(Z, full_matrices=False)
            Uz = Uz[:, :span]
            R = Nglob - Uz @ (Uz.T @ Nglob)
            sr = np.linalg.svd(R, compute_uv=False)
            missing = int(np.sum(sr > 1e-8))
            inside = float(np.abs(Bf @ Uz).max())
        else:
            span, missing, inside = 0, kdim, 0.0
        Zz = np.hstack(cols_dz)
        szz = np.linalg.svd(Zz, compute_uv=False)
        span_dz = int(np.sum(szz > 1e-10 * szz[0]))
        k = deg
        formula_patch = (2 * k - 1) * (k * nz - 1) if kind == "star" else (k - 1) * (k * nz - 1)
        formula_span = k * ncol * (k * nz - 1) if kind == "star" else (k - 1) * ncol * (k * nz - 1)
        res[kind] = {"n_patches": len(patches), "patch_size_mixed": sizes,
                     "patch_size_velocity_free": sorted(set(len(p) for p in pu)),
                     "patch_kernel_dims": sorted(set(pdims)), "formula_patch": formula_patch,
                     "span": span, "formula_span": formula_span, "missing": missing,
                     "max|B span|": inside,
                     "patch_kernel_dims_dz": sorted(set(pdims_dz)), "span_dz": span_dz}
    return res


rows = []
header = ["degree", "ncol", "nz", "L_x", "L_z", "dimV_free", "dimQ", "rank_div", "dim_ker_div", "formula_ker",
          "patch", "n_patches", "patch_size_mixed", "patch_velocity_dofs_free", "patch_kernel_dim",
          "formula_patch_kernel", "dim_span_patch_kernels", "formula_span", "kernel_dim_not_captured",
          "fraction_captured", "dim_ker_dz", "patch_kernel_dim_dz", "dim_span_patch_kernels_dz"]

configs = [
    # (nx, nz, L_x, L_z, degree)
    (2, 4, 3.0e5, 4.0e3, 2),
    (2, 6, 3.0e5, 4.0e3, 2),
    (4, 4, 3.0e5, 4.0e3, 2),
    (4, 6, 3.0e5, 4.0e3, 2),
    (4, 10, 3.0e5, 4.0e3, 2),
    (4, 6, 3.0e5, 2.5e2, 2),      # thinnest domain of the aspect-ratio sweep
    (4, 6, 3.0e5, 1.2e4, 2),      # thickest domain of the aspect-ratio sweep
    (4, 6, 1.0, 1.0, 2),          # unit square
    (2, 4, 3.0e5, 4.0e3, 1),
    (4, 6, 3.0e5, 4.0e3, 1),
    (4, 10, 3.0e5, 4.0e3, 1),
    (4, 6, 3.0e5, 2.5e2, 1),
    (4, 6, 1.0, 1.0, 1),
]

Bref = {}
for (nx, nz, length, height, deg) in configs:
    r = kernel_data(nx, nz, length, height, deg)
    report(f"degree {deg}, {r['ncol']} columns x {nz} layers, L_x = {length:g}, L_z = {height:g}: "
           f"dim V = {r['dimV']} ({r['n_bc']} constrained, {r['n_free']} free), dim Q = {r['dimQ']}, "
           f"rank(div) = {r['rank_B']} (= dim Q - 1: {r['rank_B'] == r['dimQ'] - 1}), "
           f"dim ker(div) = {r['dim_ker']} (formula k ncol (k nz - 1) + 1 = {r['formula_ker']}), "
           f"singular values at the rank gap / largest: {r['sv_gap'][0]:.2e} | {r['sv_gap'][1]:.2e}; "
           f"max|B - B_x - B_z| = {r['max|B - Bx - Bz|']:.1e}")
    for kind in ["star", "vanka"]:
        d = r[kind]
        report(f"    {kind:5s}: {d['n_patches']} patches of {d['patch_size_mixed']} DoFs of (u, u_y, b), "
               f"of which {d['patch_size_velocity_free']} free velocity DoFs; "
               f"kernel of div in a patch: dim {d['patch_kernel_dims']} (formula {d['formula_patch']}); "
               f"span of all patch kernels: dim {d['span']} (formula {d['formula_span']}), "
               f"max|B x| on the span {d['max|B span|']:.1e}; "
               f"dim of ker(div) not captured: {d['missing']} of {r['dim_ker']} "
               f"(captured fraction {(r['dim_ker'] - d['missing']) / r['dim_ker']:.4f})")
        report(f"           degraded divergence d_z u_z: dim ker = {r['dim_ker_dz']}, patch kernel dim "
               f"{d['patch_kernel_dims_dz']}, span of patch kernels {d['span_dz']}")
        rows.append([deg, r['ncol'], nz, length, height, r['n_free'], r['dimQ'], r['rank_B'], r['dim_ker'],
                     r['formula_ker'], kind, d['n_patches'], "/".join(map(str, d['patch_size_mixed'])),
                     "/".join(map(str, d['patch_size_velocity_free'])),
                     "/".join(map(str, d['patch_kernel_dims'])), d['formula_patch'], d['span'],
                     d['formula_span'], d['missing'], (r['dim_ker'] - d['missing']) / r['dim_ker'],
                     r['dim_ker_dz'], "/".join(map(str, d['patch_kernel_dims_dz'])), d['span_dz']])
    key = (nx, nz, deg)
    if key in Bref:
        B0, l0, h0 = Bref[key]
        report(f"    geometry dependence of the divergence matrix: max|B(L_x={length:g}, L_z={height:g}) - "
               f"B(L_x={l0:g}, L_z={h0:g})| = {abs(r['B'] - B0).max():.2e}   (max|B| = {abs(B0).max():.3e})")
    else:
        Bref[key] = (r["B"], length, height)

fa.write_table(os.path.join(HERE, "results_02_kernel_capture.csv"), header, rows)
with open(os.path.join(HERE, "results_02_kernel_capture.txt"), "w") as f:
    f.write("\n".join(out) + "\n")
