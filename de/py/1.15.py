from sympy.physics.units import *
from sympy import *

k = 5 *m*m/s
b = 4 *m
x, x2 = var("x, x2", positive=True)

lhs = k*6*s
tmp = b*x + x*x/2
rhs = tmp.subs(x, x2) - tmp.subs(x, 5*m)

eq = Eq( lhs, rhs )
pprint(eq)

sol = solve(eq,x2)[0]
pprint("\ns2 / m:")
tmp = sol
tmp /= m
pprint(N(tmp,3))
