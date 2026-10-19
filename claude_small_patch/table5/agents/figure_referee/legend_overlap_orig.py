'''Same overlap test for the ORIGINAL-style replicas (which reproduce the original PNGs inside the axes).'''
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import legend_overlap as LO          # patches Figure.savefig (and re-runs the restyled set first; cheap)
import replica_check as RC
LO.OUT = os.path.dirname(os.path.abspath(__file__)) + '/overlap_orig'
os.makedirs(LO.OUT, exist_ok=True)
print('---- original-style replicas ----')
for folder, name, fn in RC.jobs:
    fn(name)
