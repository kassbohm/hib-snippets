from sympy.physics.units import *
from sympy import *

b, m, g = var("b, m, g")

# 4 unknowns x'', y'', z'', lambda
x, y, z, lam = var("x, y, z, lambda ")

# 4 equations:
eq1 = Eq(m*x, lam/2)
eq2 = Eq(m*y, lam)
eq3 = Eq(m*z, -m*g + lam)
eq4 = Eq(x/2  + y + z)

sol = solve([eq1, eq2, eq3, eq4], [x, y, z, lam])

pprint(sol)

t = sqrt( 18*b/(5*g) )
x = g*t*t/9
y = 2*g*t*t/9

pprint("\nx1:")
pprint(x)
pprint("\ny1:")
pprint(y)

# ⎧   4⋅g⋅m     2⋅g     4⋅g     -5⋅g ⎫
# ⎨λ: ─────, x: ───, y: ───, z: ─────⎬
# ⎩     9        9       9        9  ⎭
#
# x1:
# 2⋅b
# ───
#  5
#
# y1:
# 4⋅b
# ───
#  5
