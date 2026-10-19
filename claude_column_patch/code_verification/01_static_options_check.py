"""
Item 1 (static part): are the PETSc options of `--patch star` exactly those of the
committed (HEAD) version of Nonlinear_Boussinesq_Irk_SC.py, and what does `--patch vanka`
change?

Method.  The solver script cannot be imported (it runs the time loop at module level),
so the statements that build the solver-parameter dictionaries are cut out of the source
with `ast` and executed in a namespace that only contains a stub of `args`.  The nested
dictionaries are then flattened the way Firedrake does it (keys joined by '_') and
compared.  No Firedrake import is needed; runs in a second.

Usage:  python 01_static_options_check.py
Output: 01_static_options_check.txt (also printed)
"""
import ast
import itertools
import os
import sys
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
COMMITTED = os.path.join(HERE, "committed_HEAD_Nonlinear_Boussinesq_Irk_SC.py")        # git show HEAD:...
WORKTREE = "/home/mh2520/molin/phd/yr2/Boussinesq_implicit/Nonlinear_Boussinesq_Irk_SC.py"

PARAM_NAMES = {"patch_params", "shifted_schur_pc_params", "params_schur"}


def names_in(node):
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def attrs_in(node):
    return {(n.value.id, n.attr) for n in ast.walk(node)
            if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name)}


def extract(path, patch, timing=False):
    """Return the flattened options of `params_schur` of the script at `path`."""
    src = open(path).read()
    tree = ast.parse(src)
    keep = []
    for stmt in tree.body:
        if isinstance(stmt, ast.Assign):
            targets = {t.id for t in stmt.targets if isinstance(t, ast.Name)}
            if targets & PARAM_NAMES:
                keep.append(stmt)
        elif isinstance(stmt, ast.If):
            test_attrs = attrs_in(stmt.test)
            touched = names_in(stmt)
            if (("args", "patch") in test_attrs or ("args", "timing") in test_attrs) and (touched & PARAM_NAMES):
                keep.append(stmt)
    module = ast.Module(body=keep, type_ignores=[])
    args = SimpleNamespace(patch=patch, timing=timing, atol=1.0e-9, rtol=1.0e-6, maxit=100, dt=50.0, shift=1.0e-8)
    ns = {"args": args, "use_rotation": False, "__name__": "__main__"}
    exec(compile(module, path, "exec"), ns)
    used = [ast.get_source_segment(src, s).splitlines()[0] for s in keep]
    return flatten(ns["params_schur"]), used


def flatten(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = prefix + str(k)
        if isinstance(v, dict):
            out.update(flatten(v, key + "_"))
        else:
            out[key] = v
    return out


def report(lines, title, a, b, name_a, name_b):
    lines.append("")
    lines.append(title)
    only_a = sorted(set(a) - set(b))
    only_b = sorted(set(b) - set(a))
    differ = sorted(k for k in set(a) & set(b) if a[k] != b[k])
    lines.append(f"  number of options: {name_a} {len(a)}, {name_b} {len(b)}")
    lines.append(f"  only in {name_a}: {len(only_a)}")
    for k in only_a:
        lines.append(f"      -{k} {a[k]!r}")
    lines.append(f"  only in {name_b}: {len(only_b)}")
    for k in only_b:
        lines.append(f"      -{k} {b[k]!r}")
    lines.append(f"  same key, different value: {len(differ)}")
    for k in differ:
        lines.append(f"      -{k}: {a[k]!r} -> {b[k]!r}")
    identical = not (only_a or only_b or differ)
    lines.append(f"  IDENTICAL: {identical}")
    return identical


lines = []
ok = True
for timing in (False, True):
    committed, used_c = extract(COMMITTED, patch=None, timing=timing)
    star, used_w = extract(WORKTREE, patch="star", timing=timing)
    vanka, _ = extract(WORKTREE, patch="vanka", timing=timing)
    if not timing:
        lines.append("statements executed from the committed script:")
        lines += ["    " + u for u in used_c]
        lines.append("statements executed from the working-tree script:")
        lines += ["    " + u for u in used_w]
    tag = "--timing" if timing else "monitored run (no --timing)"
    ok &= report(lines, f"[{tag}] committed HEAD  versus  working tree with --patch star", committed, star,
                 "committed", "star")
    report(lines, f"[{tag}] working tree --patch star  versus  --patch vanka", star, vanka, "star", "vanka")
    if not timing:
        lines.append("")
        lines.append("all options of the working tree with --patch vanka (prefix firedrake_1_ is added by Firedrake):")
        for k in sorted(vanka):
            lines.append(f"      -{k} {vanka[k]!r}")

lines.append("")
lines.append(f"RESULT: options of --patch star identical to the committed version: {ok}")
text = "\n".join(lines)
print(text)
with open(os.path.join(HERE, "01_static_options_check.txt"), "w") as f:
    f.write(text + "\n")
sys.exit(0 if ok else 1)
