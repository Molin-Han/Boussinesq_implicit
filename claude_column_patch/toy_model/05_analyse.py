'''
Analysis of the sweeps of the lowest-order toy model.

Reads   tables/spectra_full_factorial.csv   (02_spectra_sweep.py)
        tables/mg_{geom,nx,nz}_{v,full}.csv (03_mg_sweep.py)
Writes  05_analysis_summary.txt             (everything printed below)
        tables/law_test_column.csv          (condition number against the laws, per configuration)
        tables/thin_domain_test.csv         (effect of L_z at fixed N_x, N_z, dt, K)
        tables/nominal_sweeps.csv           (one-parameter sweeps around the nominal configuration)
'''
import csv
import numpy as np

NOM = dict(Nx=80, Nz=20, Lz=4000.0, dt=50.0, K=1.0)
out = open('05_analysis_summary.txt', 'w')


def log(*args):
    s = ' '.join(str(a) for a in args)
    print(s)
    out.write(s + '\n')


def read(fname):
    rows = list(csv.DictReader(open(fname)))
    for r in rows:
        for k, v in r.items():
            if k in ('patch', 'regime', 'group', 'cycle'):
                continue
            r[k] = float(v)
    return rows


def arr(rows, key):
    return np.array([r[key] for r in rows])


spec = read('tables/spectra_full_factorial.csv')
col = [r for r in spec if r['patch'] == 'column']
star = [r for r in spec if r['patch'] == 'star']
key5 = lambda r: (r['Nx'], r['Nz'], r['Lz'], r['dt'], r['K'])
spec_by = {(r['patch'],) + key5(r): r for r in spec}

log('=' * 100)
log('1. SPECTRA OF THE ONE-LEVEL PRECONDITIONED OPERATOR, FULL FACTORIAL SWEEP')
log('=' * 100)
log(f'{len(col)} configurations per patch: Nx in {sorted(set(arr(col, "Nx")))}, Nz in {sorted(set(arr(col, "Nz")))},')
log(f'Lz in {sorted(set(arr(col, "Lz")))}, dt in {sorted(set(arr(col, "dt")))}, K in {sorted(set(arr(col, "K")))}')
log(f'gamma = c/dx^2 ranges over [{arr(col, "gamma").min():.3g}, {arr(col, "gamma").max():.3g}], '
    f'eta = m dz^2/dx^2 over [{arr(col, "eta").min():.3g}, {arr(col, "eta").max():.3g}]')
for name, rows in [('star', star), ('column', col)]:
    log(f'\n{name} patch:')
    for q in ['lmax', 'lmin', 'cond', 'cond_nh', 'cond_hf', 'frac_small', 'frac_small_hf']:
        v = arr(rows, q)
        log(f'   {q:14s} min {v.min():.4g}   median {np.median(v):.4g}   max {v.max():.4g}')

# ------------------------------------------------------------------------------------------
log('\n' + '=' * 100)
log('2. LAW FOR THE CONDITION NUMBER OF THE SINGLE-COLUMN PRECONDITIONER')
log('=' * 100)
g = arr(col, 'gamma')
eta = arr(col, 'eta')
cond = arr(col, 'cond')
chf = arr(col, 'cond_hf')
lmax = arr(col, 'lmax')
lawA = g
lawB = g * np.maximum(1.5, 13.5 / (1.0 + eta))
lawC = lmax * np.maximum(g / 2.0 + 5.0 / 12.0, (4.5 * g + 1.75 + 1.375 * eta) / (1.0 + eta))
lawB_hf = 13.5 * g / (1.0 + eta)
log('law A : cond = const * gamma,                        gamma = c / dx^2')
log('law B : cond = gamma * max(1.5, 13.5 / (1 + eta)),   eta = m dz^2 / dx^2       (asymptotic, gamma >> 1)')
log('law C : cond = lmax * max(gamma/2 + 5/12, (4.5 gamma + 1.75 + 1.375 eta) / (1 + eta))   (all gamma)')
log('law B_hf : cond_hf = 13.5 gamma / (1 + eta)          (high-frequency modes only)')
log('\nrange of measured / law over the configurations with gamma above a threshold:')
log(f'{"gamma >":>10s} {"n":>6s} | {"cond/gamma":>22s} | {"cond/law B":>22s} | {"cond/law C":>22s} | {"cond_hf/law B_hf":>22s}')
for thr in [0.0, 1.0, 10.0, 100.0, 1000.0]:
    s = g > thr
    f = lambda v: f'[{v.min():8.4f}, {v.max():8.4f}]'
    log(f'{thr:10g} {int(s.sum()):6d} | {f(cond[s] / lawA[s]):>22s} | {f(cond[s] / lawB[s]):>22s} | '
        f'{f(cond[s] / lawC[s]):>22s} | {f(chf[s] / lawB_hf[s]):>22s}')

log('\nleast-squares slope of log(cond) against log(gamma) (gamma > 10), in bins of eta:')
edges = [0, 1e-3, 1e-2, 1e-1, 1, 10, 1e2, 1e3, 1e6]
for a, b in zip(edges[:-1], edges[1:]):
    s = (g > 10) & (eta >= a) & (eta < b)
    if s.sum() < 5:
        continue
    p = np.polyfit(np.log(g[s]), np.log(cond[s]), 1)
    log(f'   eta in [{a:g}, {b:g}): n = {int(s.sum()):4d}, slope {p[0]:.4f}, prefactor exp(intercept) {np.exp(p[1]):.3f}, '
        f'cond/gamma in [{(cond[s] / g[s]).min():.3f}, {(cond[s] / g[s]).max():.3f}]')

with open('tables/law_test_column.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Nx', 'Nz', 'Lz', 'dt', 'K', 'regime', 'gamma', 'eta', 'lmin', 'lmax', 'cond', 'cond_hf',
                'lawB', 'lawC', 'lawB_hf', 'cond_over_gamma', 'cond_over_lawB', 'cond_over_lawC', 'cond_hf_over_lawB_hf'])
    for i, r in enumerate(col):
        w.writerow([int(r['Nx']), int(r['Nz']), r['Lz'], r['dt'], r['K'], r['regime'], f'{g[i]:.8g}', f'{eta[i]:.8g}',
                    f'{r["lmin"]:.8g}', f'{r["lmax"]:.8g}', f'{cond[i]:.8g}', f'{chf[i]:.8g}', f'{lawB[i]:.8g}',
                    f'{lawC[i]:.8g}', f'{lawB_hf[i]:.8g}', f'{cond[i] / g[i]:.6g}', f'{cond[i] / lawB[i]:.6g}',
                    f'{cond[i] / lawC[i]:.6g}', f'{chf[i] / lawB_hf[i]:.6g}'])

log('\nin terms of the parameters of the experiments: gamma = c/dx^2 = Nx^2 max(1, 4 eps^2 m) / (4 pi^2 K)')
gg = arr(col, 'Nx') ** 2 * np.maximum(1.0, arr(col, 'four_eps2_m')) / (4 * np.pi ** 2 * arr(col, 'K'))
log(f'   max relative deviation of this formula from c/dx^2 in the table: {np.max(np.abs(gg / g - 1)):.2e}')

log('\nstar patch, for comparison:')
gs = arr(star, 'gamma')
log(f'   cond_hf in [{arr(star, "cond_hf").min():.4f}, {arr(star, "cond_hf").max():.4f}] over all {len(star)} configurations')
log(f'   eigenvalue of the uniform horizontal flow u = (1, 0): lam * (2/3 + 2 gamma) in '
    f'[{(arr(star, "lam_uniform_flow") * (2 / 3 + 2 * gs)).min():.6f}, {(arr(star, "lam_uniform_flow") * (2 / 3 + 2 * gs)).max():.6f}]'
    f'   (closed form lam = 1 / (2/3 + 2 gamma))')
log(f'   fraction of eigenvalues below 0.1 lmax: max over configurations {arr(star, "frac_small").max():.4f} (all modes), '
    f'{arr(star, "frac_small_hf").max():.4f} (high-frequency modes)')
log(f'column patch: fraction of eigenvalues below 0.1 lmax, configurations with gamma > 10: '
    f'min {arr(col, "frac_small")[g > 10].min():.4f}, max {arr(col, "frac_small")[g > 10].max():.4f} (all modes); '
    f'high-frequency modes: min {arr(col, "frac_small_hf")[g > 10].min():.4f}, max {arr(col, "frac_small_hf")[g > 10].max():.4f}')

# ------------------------------------------------------------------------------------------
log('\n' + '=' * 100)
log('3. DOES A THINNER DOMAIN (SMALLER L_z, FIXED N_x, N_z, dt, K) IMPROVE THE SINGLE-COLUMN PRECONDITIONER?')
log('=' * 100)
Lzs = sorted(set(arr(col, 'Lz')), reverse=True)
groups = sorted(set((r['Nx'], r['Nz'], r['dt'], r['K']) for r in col))
rows_out = []
stats = {}
for (Nx, Nz, dt, K) in groups:
    rr = [spec_by[('column', Nx, Nz, Lz, dt, K)] for Lz in Lzs]
    for q in ['cond', 'cond_hf']:
        v = np.array([r[q] for r in rr])
        for j in range(len(Lzs) - 1):
            # step from Lzs[j] to the thinner Lzs[j+1]
            reg = rr[j + 1]['regime'] + '->' if False else (rr[j]['regime'], rr[j + 1]['regime'])
            kk = (q, reg)
            stats.setdefault(kk, []).append(v[j + 1] / v[j])
    rows_out.append([int(Nx), int(Nz), dt, K] + [f'{r["cond"]:.6g}' for r in rr] + [f'{r["cond_hf"]:.6g}' for r in rr]
                    + [r['regime'] for r in rr])
with open('tables/thin_domain_test.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Nx', 'Nz', 'dt', 'K'] + [f'cond_Lz{Lz:g}' for Lz in Lzs] + [f'cond_hf_Lz{Lz:g}' for Lz in Lzs]
               + [f'regime_Lz{Lz:g}' for Lz in Lzs])
    w.writerows(rows_out)
log('ratio  quantity(thinner L_z) / quantity(next thicker L_z)  for every halving step of L_z (12->8 is a factor 1.5),')
log('grouped by the regimes of the two ends of the step (weak: 4 eps^2 m <= 1):')
for kk in sorted(stats):
    v = np.array(stats[kk])
    log(f'   {kk[0]:8s} {kk[1][0]:>6s} -> {kk[1][1]:<6s}: n = {v.size:5d}, ratio min {v.min():.4f}, median {np.median(v):.4f}, '
        f'max {v.max():.4f}, fraction of steps with improvement (ratio < 1): {np.mean(v < 1 - 1e-12):.3f}')
tot = {q: np.array([spec_by[('column', Nx, Nz, Lzs[-1], dt, K)][q] / spec_by[('column', Nx, Nz, Lzs[0], dt, K)][q]
                    for (Nx, Nz, dt, K) in groups]) for q in ['cond', 'cond_hf']}
for q, v in tot.items():
    log(f'   {q}: L_z = {Lzs[-1]:g} m against L_z = {Lzs[0]:g} m: ratio min {v.min():.4f}, median {np.median(v):.4f}, max {v.max():.4f}; '
        f'improved in {np.mean(v < 1 - 1e-12) * 100:.1f} % of the {v.size} (Nx, Nz, dt, K) groups')

# ------------------------------------------------------------------------------------------
log('\n' + '=' * 100)
log('4. MULTIGRID + FGMRES ITERATION COUNTS (rtol 1e-8, 3 levels, GMRES(6)/ASM smoother)')
log('=' * 100)
mg = {}
for grp in ['geom', 'nx', 'nz']:
    for cyc in ['v', 'full']:
        for r in read(f'tables/mg_{grp}_{cyc}.csv'):
            mg[(cyc, r['patch']) + key5(r)] = r


def its(cyc, patch, Nx, Nz, Lz, dt, K):
    r = mg.get((cyc, patch, float(Nx), float(Nz), float(Lz), float(dt), float(K)))
    if r is None:
        return '   -'
    return f'{int(r["its"]):4d}' if r['converged'] == 1 else '>300'


def grid(title, rowname, rowvals, colname, colvals, fixed, cyc):
    log(f'\n{title}  [{cyc} cycle]   rows: {rowname}, columns: {colname}; entries: star / column')
    log(f'   {rowname:>8s} | ' + ' | '.join(f'{colname}={c:<7g}' for c in colvals))
    for rv in rowvals:
        cells = []
        for cv in colvals:
            p = dict(fixed)
            p[rowname] = rv
            p[colname] = cv
            cells.append(its(cyc, 'star', **p).strip().rjust(4) + '/' + its(cyc, 'column', **p).strip().ljust(4))
        log(f'   {rv:8g} | ' + ' | '.join(f'{c:<{len(colname) + 8}s}' for c in cells))


DTS = sorted(set(arr(col, 'dt')))
KS = sorted(set(arr(col, 'K')))
for cyc in ['v', 'full']:
    grid('time step and shift, nominal mesh and geometry', 'dt', DTS, 'K', KS, dict(Nx=80, Nz=20, Lz=4000.0), cyc)
    grid('aspect ratio and shift, dt = 50 s', 'Lz', Lzs, 'K', KS, dict(Nx=80, Nz=20, dt=50.0), cyc)
    grid('aspect ratio and time step, K = 1', 'Lz', Lzs, 'dt', DTS, dict(Nx=80, Nz=20, K=1.0), cyc)
    grid('horizontal resolution and shift, dt = 50 s', 'Nx', [20, 40, 80, 160, 320], 'K', KS, dict(Nz=20, Lz=4000.0, dt=50.0), cyc)
    grid('vertical resolution and shift, dt = 50 s', 'Nz', [10, 20, 40], 'K', KS, dict(Nx=80, Lz=4000.0, dt=50.0), cyc)

for cyc in ['v', 'full']:
    log(f'\nsummary over all {len([k for k in mg if k[0] == cyc]) // 2} multigrid configurations, {cyc} cycle:')
    for patch in ['star', 'column']:
        rr = [r for k, r in mg.items() if k[0] == cyc and k[1] == patch]
        it = arr(rr, 'its')
        cv = arr(rr, 'converged')
        log(f'   {patch:6s}: iterations min {int(it[cv == 1].min())}, median {np.median(it):.0f}, max (converged) {int(it[cv == 1].max())}, '
            f'not converged within 300 iterations: {int(np.sum(cv == 0))}')
    for K in KS:
        rs = [r for k, r in mg.items() if k[0] == cyc and k[1] == 'star' and r['K'] == K]
        rc = [mg[(cyc, 'column') + key5(r)] for r in rs]
        a = arr(rs, 'its')
        b = arr(rc, 'its')
        log(f'   K = {K:6g}: star {int(a.min())}-{int(a.max())}, column {int(b.min())}-{int(b.max())}{"+" if np.any(arr(rc, "converged") == 0) else ""}; '
            f'column/star ratio of iterations: min {np.min(b / a):.2f}, median {np.median(b / a):.2f}, max {np.max(b / a):.2f}')

log('\nthin domain and the iteration count of the column patch (nominal mesh, geom sweep):')
for cyc in ['v', 'full']:
    ratio = []
    for dt in DTS:
        for K in KS:
            a = mg[(cyc, 'column', 80.0, 20.0, Lzs[0], dt, K)]
            b = mg[(cyc, 'column', 80.0, 20.0, Lzs[-1], dt, K)]
            ratio.append((b['its'] / a['its'], dt, K, a['its'], b['its']))
    rv = np.array([x[0] for x in ratio])
    log(f'   {cyc:4s} cycle: iterations(L_z = {Lzs[-1]:g} m) / iterations(L_z = {Lzs[0]:g} m): min {rv.min():.2f}, median {np.median(rv):.2f}, '
        f'max {rv.max():.2f}; fewer iterations on the thin domain in {int(np.sum(rv < 1))} of {rv.size} (dt, K) pairs, '
        f'equal in {int(np.sum(rv == 1))}, more in {int(np.sum(rv > 1))}')

log('\niterations of the column patch against the spectral quantity cond_hf (all multigrid configurations, converged runs):')
for cyc in ['v', 'full']:
    rr = [r for k, r in mg.items() if k[0] == cyc and k[1] == 'column' and r['converged'] == 1]
    x = np.array([spec_by[('column',) + key5(r)]['cond_hf'] for r in rr])
    y = arr(rr, 'its')
    s = x > 20
    p = np.polyfit(np.log(x[s]), np.log(y[s]), 1)
    # rank correlation
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    rho = np.corrcoef(rx, ry)[0, 1]
    log(f'   {cyc:4s} cycle: n = {x.size}, Spearman rank correlation {rho:.3f}; fit for cond_hf > 20: '
        f'iterations = {np.exp(p[1]):.3f} * cond_hf^{p[0]:.3f}')
    for lo, hi in [(0, 5), (5, 20), (20, 100), (100, 500), (500, 2500), (2500, 1e9)]:
        s = (x >= lo) & (x < hi)
        if s.sum():
            log(f'        cond_hf in [{lo:g}, {hi:g}): n = {int(s.sum()):4d}, iterations min {int(y[s].min())}, '
                f'median {np.median(y[s]):.0f}, max {int(y[s].max())}')

# nominal one-parameter sweeps in one table
with open('tables/nominal_sweeps.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['sweep', 'Nx', 'Nz', 'Lz', 'dt', 'K', 'regime', 'gamma', 'eta',
                'star_lmin', 'star_lmax', 'star_cond', 'star_cond_nh', 'star_cond_hf', 'star_n_small', 'star_n_small_hf',
                'column_lmin', 'column_lmax', 'column_cond', 'column_cond_hf', 'column_n_small', 'column_n_small_hf', 'ndof',
                'star_its_v', 'column_its_v', 'star_its_full', 'column_its_full'])
    sweeps = [('dt', 'dt', DTS), ('K', 'K', KS), ('Lz', 'Lz', Lzs), ('Nx', 'Nx', [20, 40, 80, 160, 320]), ('Nz', 'Nz', [10, 20, 40])]
    for name, key, vals in sweeps:
        for v in vals:
            p = dict(NOM)
            p[key] = v
            k5 = (float(p['Nx']), float(p['Nz']), float(p['Lz']), float(p['dt']), float(p['K']))
            s_ = spec_by[('star',) + k5]
            c_ = spec_by[('column',) + k5]
            w.writerow([name, int(p['Nx']), int(p['Nz']), p['Lz'], p['dt'], p['K'], s_['regime'], f'{s_["gamma"]:.6g}', f'{s_["eta"]:.6g}',
                        f'{s_["lmin"]:.6g}', f'{s_["lmax"]:.6g}', f'{s_["cond"]:.6g}', f'{s_["cond_nh"]:.6g}', f'{s_["cond_hf"]:.6g}',
                        int(s_['n_small']), int(s_['n_small_hf']),
                        f'{c_["lmin"]:.6g}', f'{c_["lmax"]:.6g}', f'{c_["cond"]:.6g}', f'{c_["cond_hf"]:.6g}',
                        int(c_['n_small']), int(c_['n_small_hf']), int(c_['ndof']),
                        its('v', 'star', **p).strip(), its('v', 'column', **p).strip(),
                        its('full', 'star', **p).strip(), its('full', 'column', **p).strip()])
log('\ntables/law_test_column.csv, tables/thin_domain_test.csv, tables/nominal_sweeps.csv written')
out.close()
