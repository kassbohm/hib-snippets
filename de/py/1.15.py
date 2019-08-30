# %%
from sympy.physics.units import *
from sympy import *

k = 5 *m*m/s
b = 4 *m

s2 = var("s2", positive=True)
eq = Eq( k * 6*s, b*s2 + s2*s2/2 - b *5*m - 12*m*m  )

pprint("\nEquation:")
pprint(eq)

s2 = solve(eq,s2)[0]
pprint("\ns2 / m:")
pprint(N(s2/m,3))
