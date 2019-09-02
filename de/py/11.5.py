from sympy.physics.units import *
from sympy import *

mA, mB, S, mu, g = var("mA, mB, S, mu, g")
sA, v0 = var("sA, v0")

sub_list=[
    (mA, 100*newton / (9.81 *(m/s/s) ) ),
    (mB, 200*newton / (9.81 *(m/s/s) ) ),
    (mu, 0.2),
    (g, 9.81 *(m/s/s) ),
    (v0, 1 *m/s),
    (sA, 2 *m),
]

# a = x'':
a, t = var("a, t", positive=True)

eq1 = Eq(mA * 2*a, S - mu*mA*g )
eq2 = Eq(mB * a, -2*S + mB*g )

pprint(eq1)
pprint(eq2)

sol = solve([eq1, eq2],[S,a])

pprint(sol)

aA = 2*sol[a]

# ⎧   g⋅mA⋅mB⋅(μ + 2)     -g⋅(2⋅mA⋅μ - mB) ⎫
# ⎨S: ───────────────, a: ─────────────────⎬
# ⎩      4⋅mA + mB            4⋅mA + mB    ⎭

eq = Eq(sA, aA*t*t/2 + v0*t)

sol = solve(eq, t)
t = sol[0]
t = t.simplify()
t = t.subs(sub_list)

v = aA * t + v0
v = v.subs(sub_list)
pprint("\nvA / (m/s):")
tmp = v/(m/s)
pprint(N(tmp,3))

# vA / (m/s):
# 4.68
