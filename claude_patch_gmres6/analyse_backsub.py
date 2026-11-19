'''
Figures and tables of the 2026-11-18 follow-up (make_backsub_jobs.py): the nonlinear paper sweeps with the pressure
back-substitution in the auxiliary-operator preconditioner (HEAD default since 6c4b0cf), against the 2026-11-16 sweeps
without it, and the K < 1 points rerun twice with SuperLU_DIST on the coarse level (deterministic).

    python analyse_backsub.py           # -> figures/70_backsub_nonlinear/
'''
import csv
import hashlib
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt   # noqa: E402

import analyse as A               # noqa: E402  (data access, records, plotting style; nothing runs on import)

ps, R = A.ps, A.R
CODE = 'nonlinear'
OUT = os.path.join(A.FIG, '70_backsub_nonlinear')
CONF = {   # configuration -> (label, folder function)
    'old': ('committed preconditioner (2026-11-16, MUMPS)', lambda t, p: A.folder_new(CODE, t, p)),
    'bs': ('with back-substitution (MUMPS)', lambda t, p: os.path.join(A.NEW, f'data_nonlinear_{t}_deg2_ew_gmres6_backsub_{p}')),
    'slu1': ('with back-substitution, SuperLU_DIST, run 1', lambda t, p: os.path.join(A.NEW, 'backsub_superlu', f'data_nonlinear_{t}_{p}_r1')),
    'slu2': ('with back-substitution, SuperLU_DIST, run 2', lambda t, p: os.path.join(A.NEW, 'backsub_superlu', f'data_nonlinear_{t}_{p}_r2')),
}


def digest(d, test, val, dt):
    f = A.find(d, 'residual', dt, A.key_of(test), A.target(test, val, dt))
    return hashlib.md5(open(f, 'rb').read()).hexdigest() if f else None


def collect():
    rows = []
    for t in A.TESTS:
        jt = {(c, p): A.job_table(fn(t, p)) for c, (_, fn) in CONF.items() for p in A.PATCHES}
        for val in A.members(CODE, t):
            for dt in A.SWEEP[CODE]['dts']:
                row = dict(test=t, member=val, dt=dt, **A.groups(CODE, t, val, dt))
                for c, (_, fn) in CONF.items():
                    for p in A.PATCHES:
                        if c.startswith('slu') and row['K'] >= 1:
                            continue
                        row[(c, p)] = A.record(fn(t, p), CODE, t, val, dt, jt[(c, p)])
                        if c.startswith('slu'):
                            row[('md5', c, p)] = digest(fn(t, p), t, val, dt)
                rows.append(row)
    return rows


def cnt(r):
    return r['count'] if A.ok(r) else np.nan


def kbins(rows, c):
    out = []
    for lo, hi, name in ((10, np.inf, 'K >= 10'), (1, 10, '1 <= K < 10'), (0.1, 1, '0.1 <= K < 1'), (0, 0.1, 'K < 0.1')):
        sel = [r for r in rows if lo <= r['K'] < hi and A.ok(r[(c, 'star')]) and A.ok(r[(c, 'vanka')])]
        q = [r[(c, 'vanka')]['count'] / r[(c, 'star')]['count'] for r in sel]
        e = [r[(c, 'vanka')]['count'] - r[(c, 'star')]['count'] for r in sel]
        s = [r[(c, 'star')]['count'] for r in sel]
        v = [r[(c, 'vanka')]['count'] for r in sel]
        out.append([name, str(len(sel)), f'{np.median(s):g} ({min(s)}-{max(s)})', f'{np.median(v):g} ({min(v)}-{max(v)})',
                    f'{np.median(q):.2f}', f'{max(q):.2f}', f'{np.median(e):g}', f'{max(e):g}'] if sel else [name, '0'] + [''] * 6)
    return out


def fails(rows, c, p):
    return sorted((r['test'], r['member'], r['dt']) for r in rows if (c, p) in r and r[(c, p)]['status'] not in ('ok', 'missing'))


def missing(rows):
    return sum(1 for r in rows for k, v in r.items() if isinstance(k, tuple) and len(k) == 2 and v['status'] == 'missing')


def main():
    os.makedirs(OUT, exist_ok=True)
    ps.apply()
    rows = collect()
    # ---------------------------------------------------------------- csv
    flat = []
    for r in rows:
        d = dict(test=r['test'], member=r['member'], dt=r['dt'], K=f"{r['K']:.4g}", rho1=f"{r['rho1']:.4g}")
        for c in CONF:
            for p in A.PATCHES:
                if (c, p) in r:
                    d[f'{c}_{p}'] = A.cell(r[(c, p)])
        flat.append(d)
    cols = ['test', 'member', 'dt', 'K', 'rho1'] + [f'{c}_{p}' for c in CONF for p in A.PATCHES]
    with open(os.path.join(OUT, 'counts.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, restval='')
        w.writeheader()
        w.writerows(flat)
    # ---------------------------------------------------------------- counts against dt per sweep
    for t in A.TESTS:
        mem = A.members(CODE, t)
        colors = ps.ramp(len(mem))
        dts = A.SWEEP[CODE]['dts']
        fig, axs = plt.subplots(1, 3, figsize=(14, 4.0))
        for i, val in enumerate(mem):
            rr = [r for r in rows if r['test'] == t and r['member'] == val]
            for ax, p in zip(axs[:2], A.PATCHES):
                ax.plot(dts, [cnt(r[('old', p)]) for r in rr], color=colors[i], ls=':', lw=1.0, alpha=0.7)
                ax.plot(dts, [cnt(r[('bs', p)]) for r in rr], color=colors[i], marker='o' if p == 'star' else 's', ms=3.5,
                        lw=1.4, label=A.member_tex(CODE, t, val))
            axs[2].plot(dts, [A.ratio(r[('old', 'vanka')], r[('old', 'star')]) for r in rr], color=colors[i], ls=':', lw=1.0, alpha=0.7)
            axs[2].plot(dts, [A.ratio(r[('bs', 'vanka')], r[('bs', 'star')]) for r in rr], color=colors[i], marker='o', ms=3.5, lw=1.4)
        for ax, p in zip(axs[:2], A.PATCHES):
            ax.set_xscale('log')
            ax.set_ylim(0, A.MAXIT * 1.05)
            ps.finish(ax, r'$\Delta t$ (s)', 'total FGMRES iterations (all Newton steps)', legend=(p == 'star'), small=True,
                      title=f'{A.PNAME[p]} patch')
        axs[2].set_xscale('log')
        axs[2].axhline(1.0, color=ps.GREY, lw=0.8, ls=':')
        ps.finish(axs[2], r'$\Delta t$ (s)', 'single column / star', legend=False, title='iteration ratio')
        fig.suptitle(f'nonlinear code, {t} sweep: solid = with pressure back-substitution, dotted = committed preconditioner '
                     f'of 2026-11-16 (both MUMPS, 2 ranks)', fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT, f'counts_vs_dt_{t}.png'), dpi=160)
        plt.close(fig)
    # ---------------------------------------------------------------- ratio against K
    fig, ax = plt.subplots(figsize=(7.5, 4.3))
    for c, st in (('old', dict(marker='o', mfc='none', color=ps.GREY, ms=5, ls='none', label='committed (2026-11-16)')),
                  ('bs', dict(marker='o', color='C3', ms=4, ls='none', label='with back-substitution'))):
        x = [r['K'] for r in rows if np.isfinite(A.ratio(r[(c, 'vanka')], r[(c, 'star')]))]
        y = [A.ratio(r[(c, 'vanka')], r[(c, 'star')]) for r in rows if np.isfinite(A.ratio(r[(c, 'vanka')], r[(c, 'star')]))]
        ax.plot(x, y, **st)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.axhline(1.0, color='k', lw=0.8, ls=':')
    ax.axvline(1.0, color='k', lw=0.6, ls='--')
    ps.finish(ax, r'$K = \delta/\delta_{thr}$', 'single column / star (total FGMRES)', legend=True,
              title='nonlinear code, all four sweeps (MUMPS, 2 ranks)')
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'ratio_vs_K.png'), dpi=160)
    plt.close(fig)
    # ---------------------------------------------------------------- paper-style nonlinear figures from the back-substitution runs
    for p, sub in (('star', 'star'), ('vanka', 'single_column')):
        f = lambda t: CONF['bs'][1](t, p)   # noqa: E731
        for corrected in (False, True):
            out = os.path.join(OUT, 'paper_style', sub, 'corrected' if corrected else 'paper_conventions')
            os.makedirs(out, exist_ok=True)
            for fn, *args in [
                    (R.nonlinear_history, f('dx'), 'dx', 50.0, corrected, out, 'snes_error_dx_dt50.0.png', ps.HALF),
                    (R.nonlinear_history, f('dz'), 'dz', 50.0, corrected, out, 'snes_error_dz_dt50.0.png', ps.HALF),
                    (R.nonlinear_history, f('ar'), 'ar', 100.0, corrected, out, 'snes_error_AR_dt100.0.png', ps.HALF),
                    (R.nonlinear_history, f('ar'), 'ar', 5000.0, corrected, out, 'snes_error_AR_dt5000.0.png', ps.HALF, True),
                    (R.nonlinear_history, f('dt'), 'dt', 500.0, corrected, out, 'snes_error_shift_dt500.0.png', ps.HALF),
                    (R.nonlinear_shift_sweep, f('dt'), corrected, out, 'error_shift_nonlinear.png', ps.HALF)]:
                try:
                    fn(*args)
                except Exception as e:
                    print('  skipped', type(e).__name__, e)
            plt.close('all')
    # ---------------------------------------------------------------- tables
    md = A.md_table
    hdr = ['K', 'points', 'star: median (range)', 'single column: median (range)', 'ratio median', 'ratio max',
           'extra its median', 'extra its max']
    fail_rows = []
    for c in ('old', 'bs'):
        for p in A.PATCHES:
            fl = fails(rows, c, p)
            fail_rows.append([CONF[c][0], A.PNAME[p], str(len(fl)), '; '.join(f'{t} {m:g} dt={dt:g}' for t, m, dt in fl)])
    newton = []
    for c in ('old', 'bs'):
        for p in A.PATCHES:
            n = [r[(c, p)]['newton'] for r in rows if A.ok(r[(c, p)])]
            newton.append([CONF[c][0], A.PNAME[p], f'{np.median(n):g}', f'{min(n)}-{max(n)}'] if n else [CONF[c][0], A.PNAME[p], '-', '-'])
    mesh = []
    for t, lab in (('dz', 'Nz'), ('dx', 'Nx'), ('ar', 'Lz')):
        for val in A.members(CODE, t):
            sel = [r for r in rows if r['test'] == t and r['member'] == val and r['dt'] >= 500]
            cells = []
            for c in ('old', 'bs', 'slu1'):
                for p in A.PATCHES:
                    v = [r[(c, p)]['count'] for r in sel if (c, p) in r and A.ok(r[(c, p)])]
                    cells.append(f'{np.median(v):g}' if v else '-')
            mesh.append([f'{lab} = {val:g}'] + cells)
    det, noise = [], []
    for p in A.PATCHES:
        sel = [r for r in rows if r['K'] < 1]
        same_md5 = sum(1 for r in sel if r[('md5', 'slu1', p)] and r[('md5', 'slu1', p)] == r[('md5', 'slu2', p)])
        same_cnt = sum(1 for r in sel if A.cell(r[('slu1', p)]) == A.cell(r[('slu2', p)]))
        det.append([A.PNAME[p], str(len(sel)), str(same_md5), str(same_cnt)])
        d = [r[('bs', p)]['count'] - r[('slu1', p)]['count'] for r in sel if A.ok(r[('bs', p)]) and A.ok(r[('slu1', p)])]
        noise.append([A.PNAME[p], str(len(d)), str(sum(1 for x in d if x == 0)), f'{np.median(np.abs(d)):g}' if d else '',
                      f'{max(np.abs(d)):g}' if d else ''])
    kslu = [r for r in rows if r['K'] < 1 and A.ok(r[('slu1', 'star')]) and A.ok(r[('slu1', 'vanka')])]
    q = [r[('slu1', 'vanka')]['count'] / r[('slu1', 'star')]['count'] for r in kslu]
    txt = f'''# Nonlinear sweeps with the pressure back-substitution (2026-11-18)

Generated by `claude_patch_gmres6/analyse_backsub.py`; jobs from `claude_patch_gmres6/make_backsub_jobs.py`. Paths
relative to `yr2/`.

**What was run.** The paper's four nonlinear sweeps (delta = C1/dt; 220 points per patch) for the star and the
single-column patch with `Nonlinear_Boussinesq_Irk_SC.py` at HEAD, whose auxiliary-operator preconditioner applies the
pressure back-substitution k_p <- k_p - Pi(div k_u)/delta by default since commit 6c4b0cf (2026-11-17), so that it
applies A_delta^-1 instead of T^-1 A_delta^-1 (`--no_backsub` restores the old one). The commands are byte-identical to
the 2026-11-16 sweeps (GMRES(6) smoothing, FGMRES right, Eisenstat-Walker, MUMPS coarse LU, 2 ranks); only the
script changed. Data: `data/patch_gmres6/data_nonlinear_<t>_deg2_ew_gmres6_backsub_<p>/`. In addition every point with
K < 1 (106 per patch) was run twice with SuperLU_DIST on the multigrid coarse level (deterministic):
`data/patch_gmres6/backsub_superlu/data_nonlinear_<t>_<p>_r<1,2>/`. A pilot (`data/patch_gmres6/backsub/pilot/`)
checked that HEAD's back-substitution gives the same counts as the 2026-10-16 `--pc_fix` script and that `--no_backsub`
reproduces the 2026-11-16 counts.

"committed" below = the 2026-11-16 sweeps (`../21_nonlinear_dt/` ... `../24_nonlinear_dz/`), i.e. without the
back-substitution. Counts = total FGMRES iterations of the monitored time step (all Newton steps).

## Files
- `counts_vs_dt_<t>.png`: per sweep, counts of both patches and their ratio; solid = with back-substitution, dotted = committed.
- `ratio_vs_K.png`: single column / star against K for both preconditioners.
- `counts.csv`: every point: K, rho1 and the counts (total/Newton steps) of all four configurations
  (old = committed, bs = back-substitution with MUMPS, slu1/slu2 = back-substitution with SuperLU_DIST, K < 1 only).
- `paper_style/{{star,single_column}}/{{paper_conventions,corrected}}/`: the paper's six nonlinear figures (figs. 8-12) redrawn
  from the back-substitution runs (conventions as in `../50_paper_style/README.md`).

## Failures (maxit or aborted)
Runs without output (not run or not finished): {missing(rows)}.

{md(['preconditioner', 'patch', 'number', 'points (sweep, member, dt)'], fail_rows)}

## By K (both patches converged; pooled over the four sweeps, so the nominal point counts up to 4 times)
### committed preconditioner (2026-11-16)
{md(hdr, kbins(rows, 'old'))}

### with back-substitution (MUMPS)
{md(hdr, kbins(rows, 'bs'))}

### with back-substitution, SuperLU_DIST (deterministic), K < 1 only
{md(hdr, kbins([r for r in rows if r['K'] < 1], 'slu1'))}

Single column / star with SuperLU_DIST over all {len(q)} K < 1 points: median {np.median(q) if q else np.nan:.2f}, range {min(q) if q else np.nan:.2f}-{max(q) if q else np.nan:.2f}.

## Newton steps
{md(['preconditioner', 'patch', 'median', 'range'], newton)}

## Mesh and domain dependence at large dt (medians over dt >= 500 s, where K < 1)
{md(['member', 'committed star', 'committed single col.', 'back-sub star', 'back-sub single col.', 'SuperLU star', 'SuperLU single col.'], mesh)}

## Reproducibility
SuperLU_DIST run 1 against run 2 (K < 1 points): identical residual histories (bitwise) and identical counts:
{md(['patch', 'points', 'bitwise identical', 'same count / Newton steps'], det)}

MUMPS (the back-substitution sweeps) against SuperLU_DIST run 1 at the same K < 1 points (both converged):
{md(['patch', 'points', 'same count', 'median abs. difference', 'max abs. difference'], noise)}
'''
    open(os.path.join(OUT, 'README.md'), 'w').write(txt)
    print('wrote', OUT)


if __name__ == '__main__':
    main()
