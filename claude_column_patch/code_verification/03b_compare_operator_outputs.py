"""
Compare the outputs of two runs of 03_operator_serial_vs_parallel.py DoF by DoF.

Usage: python 03b_compare_operator_outputs.py runs/operator_<tagA> runs/operator_<tagB> [more pairs ...]
       (arguments are taken two at a time)
For every stored vector y:   ||y_A - y_B||_2 / ||y_A||_2   and   max_i |y_A - y_B|_i / max_i |y_A|_i
The keys (field, base entity, x of the base entity, position in the column) of the two runs must
be identical, which is checked.
"""
import os
import sys
import numpy as np

args = sys.argv[1:]
assert len(args) >= 2 and len(args) % 2 == 0
for a, b in zip(args[0::2], args[1::2]):
    A = np.load(os.path.join(a, "outputs.npz"))
    B = np.load(os.path.join(b, "outputs.npz"))
    print("=" * 110)
    print(f"A = {a}")
    print(f"B = {b}")
    for k in sorted(A.files):
        if k.startswith('partition'):
            print(f"   {k:36s} owned base cells per rank: A {A[k].tolist()}  B {B[k].tolist()}")
    print(f"   {'vector':36s} {'n':>7s} {'||A||_2':>14s} {'rel. l2 difference':>20s} {'rel. max difference':>20s}")
    for k in sorted(A.files):
        if not k.endswith("__vals"):
            continue
        name = k[:-6]
        if k not in B.files:
            print(f"   {name:36s} missing in B")
            continue
        ka, kb = A[name + "__keys"], B[name + "__keys"]
        assert ka.shape == kb.shape and (ka == kb).all(), f"keys differ for {name}"
        ya, yb = A[k], B[k]
        d2 = np.linalg.norm(ya - yb) / np.linalg.norm(ya)
        dm = np.abs(ya - yb).max() / np.abs(ya).max()
        print(f"   {name:36s} {len(ya):7d} {np.linalg.norm(ya):14.6e} {d2:20.3e} {dm:20.3e}")
