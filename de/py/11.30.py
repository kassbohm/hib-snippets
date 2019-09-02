from sympy.physics.units import *
from sympy import *

# Given:
alpha, a, b, m, g = var("alpha, a, b, m, g", positive=True)

sub_list=[
    ( alpha,  30 *pi/180 ),
    ( a,  5 *meter ),
    ( b,  4 *meter ),
    ( g,  9.81 *meter/s**2 ),
    ]

# 4 unknowns x'', y'', z'', lambda
x, y, z, lam = var("x, y, z, lambda ")

ta, ca, sa = tan(alpha), cos(alpha), sin(alpha)

# 4 equations:
eq1 = Eq(ta*x - z)
eq2 = Eq(m*x, lam*ta)
eq3 = Eq(m*y, 0)
eq4 = Eq(m*z, -m*g - lam)

sol = solve([eq1, eq2, eq3, eq4], [x, y, z, lam])
(x, y, z, lam) = ( sol[x], sol[y], sol[z], sol[lam] )

t = var("t", positive = True)
beta, v0, xi = var("beta, v0, xi", positive=True)
tb = tan(beta)
sb, cb = sin(beta), cos(beta)

Y = var("Y")
X = x/2*xi * Y**2 + tb*ca*Y
XP = x*xi *Y + tb*ca

# 2 Equations for beta and xi:
eq1 = Eq(b*ca, X.subs(Y,a))
eq2 = Eq(0, XP.subs(Y,a))

sol = solve([eq1, eq2],[beta, xi], dict=True)[0]
pprint(sol)
xi, beta = sol[xi], sol[beta]

tmp = beta.subs(sub_list)
tmp = tmp *180/pi
pprint("\nbeta / deg:")
pprint(N(tmp,3))

cb = cos(beta)
v0 = sqrt(1/xi/cb/cb)
tmp = v0.subs(sub_list)
tmp /= meter/s
pprint("\nv0 / (m/s):")
pprint(N(tmp,3))

# ⎧       ⎛2⋅b⎞         2⋅b    ⎫
# ⎪β: atan⎜───⎟, ξ: ───────────⎪
# ⎨       ⎝ a ⎠      2         ⎬
# ⎪                 a ⋅g⋅sin(α)⎪
# ⎩                            ⎭
#
# beta / deg:
# 58.0
#
# v0 / (m/s):
# 7.39
