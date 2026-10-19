"""Move the Table-5 float block before the paragraph that discusses it (scratch copy only)."""
import sys
p = sys.argv[1]
t = open(p).read()
b = t.index("\\begin{table}[t]\n    \\centering\n    \\caption{Fieldsplit Schur-complement preconditioning versus a")
e = t.index("\\end{table}\n", b) + len("\\end{table}\n")
block = t[b:e]
t = t[:b] + t[e:]
anchor = "Table~\\ref{tab:monolithic} reports, for the geometry, time step and"
i = t.index(anchor)
t = t[:i] + block + "\n" + t[i:]
open(p, "w").write(t)
print("moved", len(block), "chars")
