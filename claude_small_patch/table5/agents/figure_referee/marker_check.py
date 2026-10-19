'''Marker-centroid comparison, original PNG vs original-style replica drawn from the restyle loader's data.
For each tab10 colour: mask, erode (removes the 1.5 pt lines, keeps the 6 pt discs), label, centroid; positions
relative to the detected axes frame. Reports, per colour, the number of markers found and the max centroid
displacement (px). Legend-box markers are excluded (only centroids inside the frame and outside the legend box).'''
import os, sys
import numpy as np
from PIL import Image
from scipy import ndimage
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import replica_check2 as RC2  # frame() and the jobs (re-renders; cheap)
RC = RC2.RC
TAB = {'C0': (31, 119, 180), 'C1': (255, 127, 14), 'C2': (44, 160, 44), 'C3': (214, 39, 40), 'C4': (148, 103, 189),
       'C5': (140, 86, 75), 'C6': (227, 119, 194)}

def centroids(img, fr, legbox):
    out = {}
    for k, c in TAB.items():
        m = np.abs(img - np.array(c)).max(axis=2) < 40
        m = ndimage.binary_erosion(m, structure=np.ones((4, 4)))
        lab, n = ndimage.label(m)
        pts = []
        for cy, cx in ndimage.center_of_mass(m, lab, range(1, n + 1)):
            x, y = cx - fr[0], cy - fr[1]
            if not (0 < x < fr[2] - fr[0] and 0 < y < fr[3] - fr[1]):
                continue
            if legbox and legbox[0] <= x <= legbox[2] and legbox[1] <= y <= legbox[3]:
                continue
            pts.append((x, y))
        out[k] = np.array(sorted(pts))
    return out

def legend_box(img, fr):
    B = img[fr[1]:fr[3], fr[0]:fr[2]]
    grey = (np.abs(B - np.array([204, 204, 204])).max(axis=2) < 12)
    ys, xs = np.nonzero(grey)
    return (xs.min() - 2, ys.min() - 2, xs.max() + 2, ys.max() + 2) if len(xs) else None

for folder, name, fn in RC.jobs:
    a = np.asarray(Image.open(f'{RC.PAPER}/{folder}/{name}').convert('RGB')).astype(int)
    b = np.asarray(Image.open(f'{RC.OUT}/{name}').convert('RGB')).astype(int)
    fa, fb = RC2.frame(a), RC2.frame(b)
    lb = legend_box(b, fb)
    ca, cb = centroids(a, fa, lb), centroids(b, fb, lb)
    parts = []
    worst = 0.0
    for k in TAB:
        A, B = ca[k], cb[k]
        if len(A) == 0 and len(B) == 0:
            continue
        if len(A) != len(B):
            parts.append(f'{k}: {len(A)} vs {len(B)} markers!'); worst = np.inf; continue
        dmax = 0.0
        for p in A:      # nearest-neighbour match
            dmax = max(dmax, np.min(np.hypot(*(B - p).T)))
        worst = max(worst, dmax)
        parts.append(f'{k}:{len(A)}')
    print(f'{name:30s} markers per colour (outside legend) {", ".join(parts)}; max centroid shift {worst:.2f} px')
