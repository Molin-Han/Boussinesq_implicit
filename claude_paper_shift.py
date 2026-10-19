'''
Shift parameter chosen from the analysis of the paper (Section "Asymptotic
shift-parameter scaling", Proposition prop:asymptotic and its dimensional form
Corollary cor:dimensional).  This module is the single source of truth for the
formula: both the command generator (Bash/bash_claude_paper.py) and the plotting
script (plot_claude_paper.py) import it, so a run and its plot can never disagree
about which delta was used.

The admissible shift is

    delta_thr = dt / (2 C_M^2),      C_M^2 = max(C_h^2, m C_v^2),

with the dimensional constants of Corollary cor:constants-dimensional,

    C_h = L_x / (2 pi),   C_v = L_z / pi,   m = lambda_max(M) = 1 + (dt N)^2 / 4 .

The two branches of the maximum are the two regimes of prop:asymptotic:

    weakly stratified   (4 eps^2 m <= 1, i.e. eps dt N <~ 1):
        delta_thr = 2 pi^2 dt / L_x^2                 ~ dt
    strongly stratified (4 eps^2 m >  1, i.e. eps dt N >~ 1):
        delta_thr = pi^2 dt / (2 m L_z^2)
                  ~ 2 pi^2 / (dt N^2 L_z^2)            ~ 1/dt   for dt N >> 2

and they meet at dt ~ 1/(eps N) = L_x / (N L_z), eps = L_z / L_x.

The experiments run with delta = K * delta_thr for a dimensionless safety factor
K.  K = 1 is the analytical (worst-case, every right-hand side) threshold, at which
every pressure mode is damped by at most 1/3.  The mode-damping factor of
eq:mode-damping for the least damped mode present in the Skamarock-Klemp data (the
(0,1) mode, Remark rem:thin-faster) is

    rho_01 = K / (K + 2 / min(1, 4 eps^2 m)),

i.e. K/(K + 2/(4 eps^2 m)) <= K/(K+2) in the weakly stratified regime and exactly
K/(K+2) in the strongly stratified one.  So at fixed K the iteration count is a
function of the single group 4 eps^2 m (~ (eps dt N)^2 once dt N >> 2) and curves
for different aspect ratios should collapse when plotted against it (and, up to the
dt N << 2 corner, against eps dt N).  That collapse is what the aspect-ratio test
checks.

The paper's *operative* threshold for these data (eq:C1-threshold-baroclinic) is
the larger delta <~ sigma_01/2; in terms of K it is K_operative = 1/min(1, 4 eps^2 m),
which is ~1400 on the nominal geometry at small dt and 1 in the strongly stratified
regime.  The K sweep of the dt test brackets it.
'''
import numpy as np

N_FREQ = 1.0e-2          # buoyancy frequency N of the test problem (1/s); = sqrt(utils.buo_freq())
LENGTH = 3.0e5           # nominal L_x (m)
HEIGHT = 4.0e3           # nominal L_z (m)


def m_factor(dt, N=N_FREQ):
    '''m = lambda_max(M) = 1 + dt^2 N^2 / 4, eq:m-J-def (M of eq:M-def).'''
    return 1.0 + (dt * N) ** 2 / 4.0


def C_h_squared(length):
    '''Horizontal Poincare constant squared, C_h = L_x / (2 pi) (periodic in x).'''
    return (length / (2.0 * np.pi)) ** 2


def C_v_squared(height):
    '''Vertical Poincare constant squared, C_v = L_z / pi (no-flux lids).'''
    return (height / np.pi) ** 2


def C_M_squared(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''C_M^2 = max(C_h^2, m C_v^2) of Corollary cor:constants-dimensional.'''
    return max(C_h_squared(length), m_factor(dt, N) * C_v_squared(height))


def regime(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''"weak" if the horizontal branch C_h^2 attains the maximum, else "strong".'''
    return 'weak' if m_factor(dt, N) * C_v_squared(height) <= C_h_squared(length) else 'strong'


def shift_threshold(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''delta_thr = dt / (2 C_M^2), eq:delta-admissible in dimensional form (s/m^2).'''
    return dt / (2.0 * C_M_squared(dt, length, height, N))


def shift_threshold_baroclinic(dt, height=HEIGHT, N=N_FREQ):
    '''
    The operative threshold for data without a depth-independent pressure component,
    sigma_01 / 2 = pi^2 dt / (2 m L_z^2) = eq:C1-threshold-baroclinic divided by dt.
    Equals shift_threshold in the strongly stratified regime.
    '''
    return dt / (2.0 * m_factor(dt, N) * C_v_squared(height))


def shift_weak_branch(dt, length=LENGTH):
    '''2 pi^2 dt / L_x^2, the weakly stratified asymptote of eq:delta-dimensional.'''
    return 2.0 * np.pi ** 2 * dt / length ** 2


def shift_strong_branch(dt, height=HEIGHT, N=N_FREQ):
    '''2 pi^2 / (dt N^2 L_z^2), the strongly stratified asymptote of eq:delta-dimensional.'''
    return 2.0 * np.pi ** 2 / (dt * N ** 2 * height ** 2)


def shift_value(K, dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''
    delta = K * delta_thr, rounded exactly as the command generator rounds it so
    that the value in the output file name can be reproduced by the plotting script.
    '''
    return float(np.round(K * shift_threshold(dt, length, height, N), decimals=20))


def group_eps_dt_N(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''The controlling group eps dt N, eps = L_z / L_x.'''
    return (height / length) * dt * N


def group_4eps2m(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''4 eps^2 m, the argument of the maximum in eq:CP-explicit; > 1 is the strong regime.'''
    return 4.0 * (height / length) ** 2 * m_factor(dt, N)


def crossover_dt(length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''dt at which the two branches meet, 1/(eps N) = L_x / (N L_z).'''
    return length / (N * height)


def K_operative(dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''
    The safety factor at which delta = K delta_thr reaches the operative threshold
    sigma_01/2 of eq:C1-threshold-baroclinic: 1/min(1, 4 eps^2 m).  Above it the (0,1)
    mode is damped by more than 1/3; the old sweeps found the GMRES boundary another
    order of magnitude or so above this.
    '''
    return 1.0 / min(1.0, group_4eps2m(dt, length, height, N))


def rho_01(K, dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''
    Predicted damping factor rho_01 = delta / (delta + sigma_01) of the (0,1)
    pressure mode at delta = K delta_thr, with sigma_01 = pi^2 dt / (m L_z^2)
    (eq:mode-eigenvalues).  In the weak regime delta_thr = dt/(2 C_h^2) and this is
    K / (K + 2/(4 eps^2 m)); in the strong regime delta_thr = sigma_01 / 2 exactly,
    so it saturates at K / (K + 2).  Together: K / (K + 2 / min(1, 4 eps^2 m)).
    '''
    return K / (K + 2.0 / min(1.0, group_4eps2m(dt, length, height, N)))


def C1_equivalent(K, dt, length=LENGTH, height=HEIGHT, N=N_FREQ):
    '''The old parametrisation delta = C_1 / dt: C_1 = delta * dt, for reference.'''
    return shift_value(K, dt, length, height, N) * dt


if __name__ == '__main__':
    print(f"Nominal geometry L_x = {LENGTH}, L_z = {HEIGHT}, eps = {HEIGHT/LENGTH}, N = {N_FREQ}, "
          f"cross-over dt = 1/(eps N) = {crossover_dt()} s")
    print(f"{'dt':>10} {'regime':>7} {'eps dt N':>10} {'4eps^2m':>10} {'delta_thr':>12} "
          f"{'weak asym':>12} {'strong asym':>12} {'C1 at K=1':>12} {'rho_01(K=1)':>12} {'K_operative':>12}")
    for dt in [1., 10., 50., 100., 500., 1e3, 5e3, 1e4, 5e4, 1e5]:
        print(f"{dt:10.0f} {regime(dt):>7} {group_eps_dt_N(dt):10.3e} {group_4eps2m(dt):10.3e} "
              f"{shift_threshold(dt):12.3e} {shift_weak_branch(dt):12.3e} {shift_strong_branch(dt):12.3e} "
              f"{C1_equivalent(1.0, dt):12.3e} {rho_01(1.0, dt):12.3e} {K_operative(dt):12.3e}")
