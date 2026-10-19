"""Referee checks of the numbers in Sec. sec:single-column against results/all_counts.csv (read only)."""
import csv, math, statistics as st
CSV = '/home/mh2520/molin/phd/yr2/Boussinesq_implicit/claude_small_patch/results/all_counts.csv'
LX, N, LZ0 = 3.0e5, 1.0e-2, 4.0e3
C1DEF = {'linear': 1e-3, 'nonlinear': 1e-5}
rows = list(csv.DictReader(open(CSV)))
D = {}
for r in rows:
    D[(r['mode'], r['test'], r['patch'], float(r['member']), float(r['dt']))] = r

def params(mode, test, member, dt):
    C1, Lz, nx, nz = C1DEF[mode], LZ0, 20, 20
    if test == 'dt': C1 = member
    if test == 'ar': Lz = member
    if test == 'dx': nx = int(member)
    if test == 'dz': nz = int(member)
    return C1, Lz, nx, nz

def K(C1, Lz, dt):
    eps = Lz / LX; m = 1 + dt**2 * N**2 / 4
    CM2 = (LX / (2 * math.pi))**2 * max(1.0, 4 * eps**2 * m)
    return (C1 / dt) / (dt / (2 * CM2))

def pairs(mode, tests=('ar', 'dt', 'dx', 'dz')):
    out = []
    for (md, test, patch, mem, dt), r in D.items():
        if md != mode or patch != 'star' or test not in tests: continue
        v = D[(md, test, 'vanka', mem, dt)]
        C1, Lz, nx, nz = params(md, test, mem, dt)
        out.append(dict(test=test, mem=mem, dt=dt, C1=C1, Lz=Lz, nx=nx, nz=nz, K=K(C1, Lz, dt),
                        s=r['status'], v=v['status'],
                        cs=float(r['count']) if r['status'] == 'ok' else None,
                        cv=float(v['count']) if v['status'] == 'ok' else None,
                        ns=r['newton_steps'], nv=v['newton_steps']))
    return out

def summ(xs):
    if not xs: return 'none'
    return f'n={len(xs)} median {st.median(xs):.2f} min {min(xs):.2f} max {max(xs):.2f}'

for mode in ('linear', 'nonlinear'):
    P = pairs(mode)
    both = [p for p in P if p['cs'] is not None and p['cv'] is not None]
    print(f'\n===== {mode}: {len(P)} pairs, {len(both)} both converged')
    print('non-converged star:', sorted((p['test'], p['mem'], p['dt'], p['s']) for p in P if p['cs'] is None))
    print('non-converged col :', sorted((p['test'], p['mem'], p['dt'], p['v']) for p in P if p['cv'] is None))
    print('star fails, column converges:', [(p['test'], p['mem'], p['dt'], p['cv']) for p in P if p['cs'] is None and p['cv'] is not None])
    print('column fails, star converges:', [(p['test'], p['mem'], p['dt'], p['cs']) for p in P if p['cv'] is None and p['cs'] is not None])
    if mode == 'linear':
        big = [p for p in both if abs(p['cv'] - p['cs']) > 1]
        print(f'|col-star|>1 in {len(big)} of {len(both)} runs:')
        for p in sorted(big, key=lambda p: (p['test'], p['mem'], p['dt'])):
            print(f"  {p['test']:2s} mem={p['mem']:<8g} dt={p['dt']:<7g} star={p['cs']:.0f} col={p['cv']:.0f} ratio={p['cv']/p['cs']:.2f} K={p['K']:.3g}")
        print('within one iteration:', sum(abs(p['cv'] - p['cs']) <= 1 for p in both), 'of', len(both))
        print('ratio all:', summ([p['cv'] / p['cs'] for p in both]))
        nonthin = [p for p in both if p['Lz'] >= 4000]
        print('ratio L_z >= 4 km (incl dt sweep):', summ([p['cv'] / p['cs'] for p in nonthin]))
        print('max ratio by test:', {t: max(p['cv'] / p['cs'] for p in both if p['test'] == t) for t in ('ar', 'dt', 'dx', 'dz')})
    else:
        for lab, sel in (('K>=1', lambda p: p['K'] >= 1), ('K<1', lambda p: p['K'] < 1)):
            Q = [p for p in both if sel(p)]
            print(f'{lab}: ratio {summ([p["cv"]/p["cs"] for p in Q])}; identical {sum(p["cv"]==p["cs"] for p in Q)} of {len(Q)}')
            for t in ('ar', 'dt', 'dx', 'dz'):
                Qt = [p for p in Q if p['test'] == t]
                print(f'    {t}: {summ([p["cv"]/p["cs"] for p in Qt])}; identical {sum(p["cv"]==p["cs"] for p in Qt)} of {len(Qt)}')
        # C1 = 1e-5 runs only (ar, dx, dz + dt column)
        c15 = [p for p in both if abs(p['C1'] - 1e-5) < 1e-12]
        for lab, sel in (('dt<=10', lambda p: p['dt'] <= 10), ('dt 50-100', lambda p: 50 <= p['dt'] <= 100), ('dt>=500', lambda p: p['dt'] >= 500), ('dt<=100', lambda p: p['dt'] <= 100)):
            Q = [p for p in c15 if sel(p)]
            print(f'C1=1e-5 {lab}: ratio {summ([p["cv"]/p["cs"] for p in Q])}; identical {sum(p["cv"]==p["cs"] for p in Q)} of {len(Q)}; star {min(p["cs"] for p in Q):.0f}-{max(p["cs"] for p in Q):.0f}, col {min(p["cv"] for p in Q):.0f}-{max(p["cv"] for p in Q):.0f}')
        print('C1=1e-5, dt in {50,100}: per run')
        for p in sorted([p for p in c15 if p['dt'] in (50.0, 100.0)], key=lambda p: (p['dt'], p['test'], p['mem'])):
            print(f"  dt={p['dt']:<5g} {p['test']} mem={p['mem']:<7g} star={p['cs']:.0f}/{p['ns']} col={p['cv']:.0f}/{p['nv']} ratio={p['cv']/p['cs']:.2f} K={p['K']:.2f}")
        print('C1=1e-5, dt<=10 not identical:')
        for p in [p for p in c15 if p['dt'] <= 10 and p['cv'] != p['cs']]:
            print(f"  dt={p['dt']:<5g} {p['test']} mem={p['mem']:<7g} star={p['cs']:.0f}/{p['ns']} col={p['cv']:.0f}/{p['nv']}")
        print('K>=1 runs with ratio >= 1.3:')
        for p in sorted([p for p in both if p['K'] >= 1 and p['cv'] / p['cs'] >= 1.3], key=lambda p: p['cv']/p['cs']):
            print(f"  {p['test']} mem={p['mem']:<7g} dt={p['dt']:<6g} C1={p['C1']:.0e} star={p['cs']:.0f} col={p['cv']:.0f} ratio={p['cv']/p['cs']:.2f} K={p['K']:.2f}")
        print('K<1 runs with ratio <= 1.45:')
        for p in sorted([p for p in both if p['K'] < 1 and p['cv'] / p['cs'] <= 1.45], key=lambda p: p['cv']/p['cs']):
            print(f"  {p['test']} mem={p['mem']:<7g} dt={p['dt']:<6g} C1={p['C1']:.0e} star={p['cs']:.0f} col={p['cv']:.0f} ratio={p['cv']/p['cs']:.2f} K={p['K']:.3f}")
        print('max column count, ar/dx/dz, dt>=50:', max(p['cv'] for p in both if p['test'] != 'dt' and p['dt'] >= 50))
        for nz in (10, 20, 40, 60, 80):
            Q = [p for p in both if p['test'] == 'dz' and p['nz'] == nz and p['dt'] >= 500]
            print(f'  dz N_z={nz} dt>=500: col {[int(p["cv"]) for p in sorted(Q, key=lambda p: p["dt"])]} star {[int(p["cs"]) for p in sorted(Q, key=lambda p: p["dt"])]}')
        for nx in (10, 20, 40, 60, 80):
            Q = [p for p in both if p['test'] == 'dx' and p['nx'] == nx and p['dt'] >= 500]
            print(f'  dx N_x={nx} dt>=500: col {[int(p["cv"]) for p in sorted(Q, key=lambda p: p["dt"])]} star {[int(p["cs"]) for p in sorted(Q, key=lambda p: p["dt"])]}')
        # faster where ratio <= 1/0.7
        r = [p['cv'] / p['cs'] for p in both]
        print('ratio <= 1.43 in', sum(x <= 1/0.7 for x in r), 'of', len(r), '; K<1 & ratio<=1.43:', sum(1 for p in both if p['K'] < 1 and p['cv']/p['cs'] <= 1/0.7), 'of', sum(1 for p in both if p['K'] < 1))

print('\n===== nonlinear ratio by K bin (all sweeps, both converged)')
P = [p for p in pairs('nonlinear') if p['cs'] is not None and p['cv'] is not None]
for lo, hi in ((100, 1e99), (10, 100), (1, 10), (0.1, 1), (0, 0.1)):
    Q = [p for p in P if lo <= p['K'] < hi]
    rr = [p['cv'] / p['cs'] for p in Q]
    print(f'K in [{lo:g},{hi:g}): {summ(rr)}; identical {sum(p["cv"]==p["cs"] for p in Q)} of {len(Q)}; dts {sorted(set(p["dt"] for p in Q))}')
print('\n===== linear ratio by K bin')
P = [p for p in pairs('linear') if p['cs'] is not None and p['cv'] is not None]
for lo, hi in ((100, 1e99), (10, 100), (1, 10), (0.1, 1), (0, 0.1)):
    Q = [p for p in P if lo <= p['K'] < hi]
    rr = [p['cv'] / p['cs'] for p in Q]
    print(f'K in [{lo:g},{hi:g}): {summ(rr)}; within one {sum(abs(p["cv"]-p["cs"])<=1 for p in Q)} of {len(Q)}')
print('\nnonlinear ar sweep, column counts by L_z at dt>=500:')
P = [p for p in pairs('nonlinear', ('ar',)) if p['cv'] is not None]
for Lz in (12000, 8000, 4000, 2000, 1000):
    Q = sorted([p for p in P if p['Lz'] == Lz and p['dt'] >= 500], key=lambda p: p['dt'])
    print(f'  L_z={Lz}: col {[int(p["cv"]) for p in Q]} star {[int(p["cs"]) for p in Q]} median col {st.median([p["cv"] for p in Q])}')
