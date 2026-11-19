'''
Counts of the nonlinear sweeps in the paper's shift convention (make_paperconv_jobs.py; --shift C1/(2 dt), so that
delta = C1/dt as in the linear runs), both patches, and the numbers the paper quotes from its own (Mac) runs of the
same commands, recomputed from these runs.

    MPLBACKEND=Agg python analyse_paperconv.py      -> figures/80_paperconv_nonlinear/{README.md,counts.csv}
'''
import csv
import math
import os

import numpy as np

import analyse as A

CODE = 'nonlinear'
OUT = os.path.join(A.FIG, '80_paperconv_nonlinear')
DATA = os.path.join(A.NEW, 'paperconv')
LX = 3.0e5
S = A.SWEEP[CODE]


def rec(test, p, val, dt):
    d = os.path.join(DATA, f'data_nonlinear_{test}_{p}')
    key, tgt = A.key_of(test), A.target(test, val, dt)
    r = dict(count=np.nan, newton=np.nan, status='missing')
    cf = A.find(d, 'snes_ksp_cum', dt, key, tgt / 2 if key == 'shift' else tgt)
    if cf is not None:
        cum = np.atleast_1d(np.loadtxt(cf)).astype(int)
        per = np.diff(cum) if len(cum) > 1 else cum
        r.update(count=int(cum[-1]), newton=len(cum) - 1, status='maxit' if np.any(per >= A.MAXIT) else 'ok')
    return r


def c1_of(test, val, dt):
    return val if test == 'dt' else S['C1']


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for t in A.TESTS:
        for val in A.members(CODE, t):
            for dt in S['dts']:
                c1 = c1_of(t, val, dt)
                rows.append(dict(test=t, member=val, dt=dt, C1=c1, admissible=c1 <= 2 * math.pi ** 2 * dt ** 2 / LX ** 2,
                                 star=rec(t, 'star', val, dt), vanka=rec(t, 'vanka', val, dt)))
    with open(os.path.join(OUT, 'counts.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['test', 'member', 'dt', 'C1_paper', 'admissible', 'star_count', 'star_newton', 'star_status',
                    'vanka_count', 'vanka_newton', 'vanka_status'])
        for r in rows:
            w.writerow([r['test'], r['member'], r['dt'], r['C1'], int(r['admissible'])] +
                       [r[p][k] for p in ('star', 'vanka') for k in ('count', 'newton', 'status')])
    L = []
    pr = L.append
    ok = lambda r: r['star']['status'] == 'ok' and r['vanka']['status'] == 'ok'   # noqa: E731
    miss = [(r['test'], r['member'], r['dt'], p) for r in rows for p in ('star', 'vanka') if r[p]['status'] == 'missing']
    pr(f'runs: {len(rows)} per patch; missing outputs: {len(miss)} {miss[:6]}')
    for p in ('star', 'vanka'):
        f = [(r['test'], r['member'], r['dt']) for r in rows if r[p]['status'] == 'maxit']
        pr(f'{p} failures ({len(f)}): {f}')
    sel = lambda t, dt: [r for r in rows if r['test'] == t and r['dt'] == dt]   # noqa: E731
    for t in ('dx', 'dz'):
        pr(f"{t} dt=50: " + ', '.join(f"{r['member']}: {r['star']['count']}/{r['star']['newton']} col {r['vanka']['count']}/{r['vanka']['newton']}" for r in sel(t, 50.)))
    pr('ar dt=100: ' + ', '.join(f"{r['member']}: {r['star']['count']}/{r['star']['newton']} col {r['vanka']['count']}/{r['vanka']['newton']}" for r in sel('ar', 100.)))
    pr('dt sweep (star count/newton, col count) per C1:')
    for c in S['C1_list']:
        pr(f'  C1={c:g}: ' + ', '.join(f"{r['dt']:g}:{r['star']['count'] if r['star']['status']=='ok' else r['star']['status']}/{r['star']['newton']} ({r['vanka']['count'] if r['vanka']['status']=='ok' else r['vanka']['status']})"
                                        for r in rows if r['test'] == 'dt' and r['member'] == c))
    # lowest C1 check
    for dt in S['dts']:
        d = {r['member']: r['star']['count'] for r in rows if r['test'] == 'dt' and r['dt'] == dt and r['star']['status'] == 'ok'}
        if not d:
            continue
        lo = min(d.values())
        pr(f'  dt={dt:g}: min star {lo} at C1 {[c for c in d if d[c] == lo]}; C1<=1e-5: {[d.get(c) for c in (1e-5, 5e-6, 1e-6)]}')
    nw = [r[p]['newton'] for r in rows for p in ('star', 'vanka') if r[p]['status'] == 'ok']
    pr(f'Newton steps over converged runs: {sorted(set(nw))}; star only: {sorted(set(r["star"]["newton"] for r in rows if r["star"]["status"] == "ok"))}')
    both = [r for r in rows if ok(r)]
    adm = [r for r in both if r['admissible']]
    oth = [r for r in both if not r['admissible']]
    q = lambda rs: [r['vanka']['count'] / r['star']['count'] for r in rs]   # noqa: E731
    pr(f'both converged: {len(both)}; admissible: {len(adm)}, ratio median {np.median(q(adm)):.2f}, max {max(q(adm)):.2f}, '
       f'col median {np.median([r["vanka"]["count"] for r in adm]):g}, star median {np.median([r["star"]["count"] for r in adm]):g}')
    pr(f'  other: {len(oth)}, ratio median {np.median(q(oth)):.2f}, equal in {sum(r["vanka"]["count"] == r["star"]["count"] for r in oth)}, '
       f'max {max(q(oth)):.2f}')
    for n in S['nzs']:
        a = [r for r in adm if r['test'] == 'dz' and r['member'] == n]
        pr(f'  Nz={n} admissible: star median {np.median([r["star"]["count"] for r in a]):g}, col median {np.median([r["vanka"]["count"] for r in a]):g} ({len(a)} runs)')
    big = sorted(adm, key=lambda r: -r['vanka']['count'] / r['star']['count'])[:5]
    pr('  largest admissible ratios: ' + '; '.join(f"{r['test']} {r['member']} dt={r['dt']:g}: {r['vanka']['count']}/{r['star']['count']}" for r in big))
    txt = '\n'.join(L)
    print(txt)
    open(os.path.join(OUT, 'checks.txt'), 'w').write(txt + '\n')


if __name__ == '__main__':
    main()
