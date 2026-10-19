'''Quantitative version of replica_check.py: compare original PNG and original-style replica (drawn from the data that
restyle_paper_figures.py plots) INSIDE THE AXES FRAME ONLY, with the legend box masked out, after locating the axes
frame in each image (so tight-bbox crops and font differences do not matter). Zero differing pixels = same curves,
same end points, same axis ranges.'''
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replica_check as RC   # re-uses the drawing functions (this re-runs replica_check once; cheap)

OUT = RC.OUT
PAPER = RC.PAPER

def frame(img):
    '''axes frame from the long black spine runs: returns (x0, y0, x1, y1) of the spines.'''
    g = img.mean(axis=2)
    blk = g < 80
    rows = np.nonzero(blk.sum(axis=1) > 0.5 * img.shape[1])[0]
    cols = np.nonzero(blk.sum(axis=0) > 0.5 * img.shape[0])[0]
    return cols.min(), rows.min(), cols.max(), rows.max()

def legend_box_rel(name):
    '''legend extent relative to the axes frame, from a fresh replica draw (pixels at dpi 100).'''
    return None

results = []
for folder, name, fn in RC.jobs:
    a = np.asarray(Image.open(f'{PAPER}/{folder}/{name}').convert('RGB')).astype(int)
    b = np.asarray(Image.open(f'{OUT}/{name}').convert('RGB')).astype(int)
    fa, fb = frame(a), frame(b)
    wa, ha = fa[2] - fa[0], fa[3] - fa[1]
    wb, hb = fb[2] - fb[0], fb[3] - fb[1]
    if (wa, ha) != (wb, hb):
        print(f'{name:32s} FRAME SIZE DIFFERS {wa}x{ha} vs {wb}x{hb}'); continue
    m = 3   # skip spines and inward tick marks
    A = a[fa[1] + m:fa[3] - m, fa[0] + m:fa[2] - m]
    B = b[fb[1] + m:fb[3] - m, fb[0] + m:fb[2] - m]
    d = np.abs(A - B).max(axis=2) > 40
    # legend: the region where the original has its light-grey legend frame; find it as the bounding box of the
    # differing pixels' cluster that coincides with the legend in the replica (white box with grey edge, alpha 0.8)
    # simpler: mask the rectangle spanned by the legend frame colour (204,204,204)-ish in the replica
    grey = (np.abs(B - np.array([204, 204, 204])).max(axis=2) < 12)
    ys, xs = np.nonzero(grey)
    mask = np.zeros_like(d)
    if len(xs):
        mask[ys.min():ys.max() + 1, xs.min():xs.max() + 1] = True
    outside = d & ~mask
    yy, xx = np.nonzero(outside)
    print(f'{name:32s} frame {wa}x{ha}px; differing pixels inside axes: {d.sum():5d}; outside legend box: {outside.sum():4d}'
          + (f' at x {xx.min()}-{xx.max()}, y {yy.min()}-{yy.max()}' if outside.any() else ''))
    Image.fromarray(np.where(outside[..., None], np.array([255, 0, 0]), A).astype(np.uint8)).save(f'{OUT}/axes_outside_legend_diff_{name}')
