"""
01_kernel_dimension_counts.py

Numerical check of the dimension counts of Sections 2-4 of theory_kernel.md.
For k = 1, 2 and a few small meshes it computes (by SVD rank computations on the
tensor-product model in fe1d.py)

  dim V_h, dim N_h = dim ker(div), dim curl(Psi_h), whether (1,0) is in N_h
  for the star and the single-column (vanka) patch:
      dim V_j, dim (N_h ∩ V_j), dim sum_j (N_h ∩ V_j),
      dim [ N_H + sum_j (N_h ∩ V_j) ]     (coarse level: coarsening by 2 in x)
  and compares every number with the closed formula.

Run:  python 01_kernel_dimension_counts.py  > 01_kernel_dimension_counts.txt
"""
import numpy as np
from fe1d import Model2D, rank, nullspace


def kernel_in_subspace(B, idx, n):
    """basis (as columns of length n) of { u in span(e_idx) : B u = 0 }"""
    Bj = B[:, idx]
    Nj = nullspace(Bj)
    out = np.zeros((n, Nj.shape[1]))
    out[idx, :] = Nj
    return out


def main():
    print("Dimension counts: computed (SVD rank) versus formula")
    print("=" * 100)
    allok = True
    for k in (1, 2):
        for (Nx, nz) in ((8, 3), (8, 5), (12, 4)):
            Lx, Lz = 3.0, 1.0
            fine = Model2D(Nx, nz, Lx, Lz, k)
            coarse = Model2D(Nx // 2, nz, Lx, Lz, k)
            B = fine.div_matrix()
            n = fine.n
            # --- global kernel
            N = nullspace(B)
            dimN = N.shape[1]
            f_dimV = k * Nx * (2 * k * nz - 1)
            f_dimN = k * Nx * (k * nz - 1) + 1
            C = fine.curl_matrix()
            r_curl = rank(C)
            div_of_curl = np.abs(B @ C).max()
            e = fine.harmonic()
            div_e = np.abs(B @ e).max()
            e_in_curl = rank(np.hstack([C, e[:, None]])) == r_curl
            print(f"\nk = {k}, Nx = {Nx}, nz = {nz}")
            print(f"  dim V_h              = {n:5d}   formula k Nx (2 k nz - 1)        = {f_dimV}")
            print(f"  rank div             = {rank(B):5d}   formula k^2 Nx nz - 1            = {k*k*Nx*nz - 1}")
            print(f"  dim N_h              = {dimN:5d}   formula k Nx (k nz - 1) + 1      = {f_dimN}")
            print(f"  dim curl(Psi_h)      = {r_curl:5d}   formula k Nx (k nz - 1)          = {k*Nx*(k*nz-1)}"
                  f"   (max |div curl| = {div_of_curl:.1e})")
            print(f"  (1,0): max|div| = {div_e:.1e}, contained in curl(Psi_h): {e_in_curl}")
            allok &= (n == f_dimV) and (dimN == f_dimN) and (r_curl == k * Nx * (k * nz - 1)) and (not e_in_curl)

            # coarse kernel, prolonged
            P = fine.prolongation_from(coarse)
            NH = P @ nullspace(coarse.div_matrix())
            assert np.abs(B @ NH).max() < 1e-10
            print(f"  dim N_H (coarse)     = {NH.shape[1]:5d}   formula k (Nx/2) (k nz - 1) + 1  = {k*(Nx//2)*(k*nz-1)+1}")

            for kind, f_dimVj, f_dimNj, f_sum, f_sumH in (
                ("star",
                 k * (4 * k - 1) * nz - 2 * k,
                 (2 * k - 1) * (k * nz - 1),
                 k * Nx * (k * nz - 1),
                 k * Nx * (k * nz - 1) + 1),
                ("vanka",
                 k * (2 * k + 1) * nz - k,
                 (k - 1) * (k * nz - 1),
                 (k - 1) * Nx * (k * nz - 1),
                 # coarse kernel + fine bubbles: k = 1 -> only N_H ; k >= 2 -> everything
                 (k * (Nx // 2) * (k * nz - 1) + 1) if k == 1 else (k * Nx * (k * nz - 1) + 1)),
            ):
                pats = fine.patches(kind)
                dimVj = len(pats[0])
                Kj = [kernel_in_subspace(B, idx, n) for idx in pats]
                dimNj = Kj[0].shape[1]
                allK = np.hstack(Kj)
                s = rank(allK) if allK.shape[1] > 0 else 0
                sH = rank(np.hstack([NH, allK]))
                print(f"  {kind:5s}: dim V_j = {dimVj:4d} (formula {f_dimVj:4d}) | "
                      f"dim N_h∩V_j = {dimNj:3d} (formula {f_dimNj:3d}) | "
                      f"dim sum_j N_h∩V_j = {s:4d} (formula {f_sum:4d}) | "
                      f"dim N_H + sum_j = {sH:4d} (formula {f_sumH:4d}) | dim N_h = {dimN}")
                allok &= (dimVj == f_dimVj) and (dimNj == f_dimNj) and (s == f_sum) and (sH == f_sumH)
    print("\nAll counts agree with the formulas:", allok)

    print("\nNominal configuration (Nx = 80, nz = 20, k = 2), formulas only:")
    k, Nx, nz = 2, 80, 20
    print(f"  dim V_h = {k*Nx*(2*k*nz-1)}, dim N_h = {k*Nx*(k*nz-1)+1}")
    print(f"  star : dim V_j = {k*(4*k-1)*nz-2*k}, dim N_h∩V_j = {(2*k-1)*(k*nz-1)}, sum = {k*Nx*(k*nz-1)}")
    print(f"  vanka: dim V_j = {k*(2*k+1)*nz-k}, dim N_h∩V_j = {(k-1)*(k*nz-1)}, sum = {(k-1)*Nx*(k*nz-1)}")


if __name__ == "__main__":
    main()
