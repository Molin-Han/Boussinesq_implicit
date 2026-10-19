from firedrake import Constant
c = Constant(0.0)
try:
    print('bool(Constant(0.0)) =', bool(c))
except Exception as e:
    print('bool raises', type(e).__name__, e)
