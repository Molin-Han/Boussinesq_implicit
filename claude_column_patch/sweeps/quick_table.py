'''
Quick look at a (possibly unfinished) sweep folder: total KSP iterations (last entry of
snes_ksp_cum_*.out) and the per-Newton-step counts, read directly from the .out files.

    python quick_table.py <data folder> <test: dx|dz|ar|dt>
'''
import sys
import glob
import os
import re
import numpy as np

folder, test = sys.argv[1], sys.argv[2]
key = {'dx': 'nx', 'dz': 'nz', 'ar': 'ar', 'dt': 'shift'}[test]
data = {}
for f in glob.glob(os.path.join(folder, f'snes_ksp_cum_dt*_{key}*.out')):
    m = re.match(rf'snes_ksp_cum_dt([0-9.e+-]+)_{key}([0-9.e+-]+)\.out', os.path.basename(f))
    dt, val = float(m.group(1)), float(m.group(2))
    cum = np.atleast_1d(np.loadtxt(f))
    data[(dt, val)] = cum
dts = sorted({k[0] for k in data})
vals = sorted({k[1] for k in data})
if test == 'dt':
    # ! the shift differs for every dt: print one line per dt with the shifts in increasing order
    for dt in dts:
        row = sorted((v, data[(dt, v)]) for (d, v) in data if d == dt)
        print(f'dt={dt:g}: ' + '  '.join(f'{v:.3g}:{int(c[-1])}' for v, c in row))
else:
    print(f'{"dt":>9} | ' + ' | '.join(f'{key}={v:g}'.rjust(14) for v in vals))
    for dt in dts:
        cells = []
        for v in vals:
            c = data.get((dt, v))
            cells.append('-'.rjust(14) if c is None else (f'{int(c[-1])} (' + '/'.join(str(int(x)) for x in np.diff(c)) + ')').rjust(14))
        print(f'{dt:9g} | ' + ' | '.join(cells))
