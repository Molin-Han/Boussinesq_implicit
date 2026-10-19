# First-kind Chebyshev error polynomial of degree k on [a,b] (PETSc KSPCHEBYSHEV, fixed interval)
import numpy as np
from numpy.polynomial import chebyshev as C
a, b, k = 0.2, 2.2, 6
Tk = lambda s: C.chebval(s, [0]*k + [1])
q = lambda t: Tk((b + a - 2*t)/(b - a)) / Tk((b + a)/(b - a))
t = np.linspace(a, b, 20001)
print(f'max|q| on [a,b] = {abs(q(t)).max():.4e}  (=1/T_6(1.2) = {1/Tk(1.2):.4e})')
print(f'q(2) = {q(2.0):.4e}  q(2)^2 = {q(2.0)**2:.4e}')
h = 1e-7; print(f"q'(0) = {(q(h)-q(-h))/(2*h):.4f}")
mu = np.array([4.3847e-04, 1.7524e-03, 3.9421e-03, 7.0215e-03])
pred = 1 - (q(2*mu)**2) * (1 - mu)**2
print('barotropic mu:', mu, '\n predicted MONO eig 1-q(2mu)^2(1-mu)^2:', np.round(pred, 4), '\n ratio pred/mu:', np.round(pred/mu, 1))
for m in [0.33, 0.66, 0.82, 0.9, 0.976]:
    print(f'baroclinic mu={m}: B_2A~2mu={2*m:.2f} -> |q(2mu)|^2 = {q(2*m)**2:.2e}')
