'''Candidate replacements (written to ./cand/ only) for the restyled figures whose legend hides data, plus the
optional true-count / converged-point variants. Uses restyle_paper_figures.py unchanged; the legend is rebuilt just
before saving (see legend_fix_trial.py), and each candidate is run through the same overlap test.'''
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import legend_fix_trial as T          # defines patched savefig, legend(), ylim_top(), relabel_short_ar(), chain()
R, ps = T.R, T.ps
T.OUT = os.path.dirname(os.path.abspath(__file__)) + '/cand'
os.makedirs(T.OUT, exist_ok=True)
nolegend = lambda fig: fig.axes[0].get_legend().remove()
_rep = T.overlap_report
T.overlap_report = lambda fig, name: 'no legend' if fig.axes[0].get_legend() is None else _rep(fig, name)
keep = lambda fig: None
T.TRIAL.clear()
T.TRIAL['error_AR_Robust_dt50.0.png'] = [('shortlabels', T.chain(T.relabel_short_ar, T.legend(loc='upper right')))]
T.TRIAL['error_shift.png'] = [('ncol3_ylim120', T.chain(T.legend(loc='upper right', ncol=3, columnspacing=1.0), T.ylim_top(120)))]
T.TRIAL['snes_error_shift_dt500.0.png'] = [('nolegend', nolegend)]
T.TRIAL['snes_error_AR_dt5000.0.png'] = [('keepconverged', keep)]
ps.apply()
d = R.folder('paper', 'linear', 'ar');  R.linear_history(d, 'ar', 50.0, False, T.OUT, 'error_AR_Robust_dt50.0.png', ps.SINGLE)
d = R.folder('paper', 'linear', 'dt');  R.linear_shift_sweep(d, False, T.OUT, 'error_shift.png', ps.SINGLE)
T.TRIAL['error_shift.png'] = [('corrected_ncol3_ylim120', T.chain(T.legend(loc='upper right', ncol=3, columnspacing=1.0), T.ylim_top(120)))]
R.linear_shift_sweep(d, True, T.OUT, 'error_shift.png', ps.SINGLE)
d = R.folder('paper', 'nonlinear', 'dt'); R.nonlinear_history(d, 'dt', 500.0, False, T.OUT, 'snes_error_shift_dt500.0.png', ps.HALF)
d = R.folder('paper', 'nonlinear', 'ar'); R.nonlinear_history(d, 'ar', 5000.0, False, T.OUT, 'snes_error_AR_dt5000.0.png', ps.HALF, False)
