from sympy.physics.units import *
from sympy import *

prec = 3
g = var("g")
cm = m/100
newton = kg*m/s**2

R, M, I, F, t = var("R, M, I, F, t ")
S, T, x, p = var("S, T, x, p")

sub_list=[
    (R, 10 *cm),
    (M, 2 *kg),
    (I, 2*kg*m**2 / 10),
    (F, 40 *newton),
    (t, 2*s),
    (g, 981 *m/s**2 / 100),
    ]

eq1 = Eq( 2*I/R * x, R*(F - S) )
eq2 = Eq(  M * x, 2*S - M*g )

sol = solve( [eq1, eq2], [S, x])

pprint("\nx'':")
x = sol[x]
pprint(x)
pprint("\nx'' / (m/s²):")
tmp = x
tmp = tmp.subs(sub_list)
tmp /= (m/s**2)
pprint(N(tmp,prec))

pprint("\nv / (m/s):")
v = x*t
tmp = v
tmp = tmp.subs(sub_list)
tmp /= m/s
pprint(N(tmp,prec))

# x'':
#  2
# R ⋅(2⋅F - M⋅g)
# ──────────────
#            2
#   4⋅I + M⋅R
#
# x'' / (m/s²):
# 0.736
#
# v / (m/s):
# 1.47
