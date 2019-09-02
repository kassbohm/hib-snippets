from sympy.physics.mechanics import *
from sympy.physics.units import *
from sympy import *

g, m, c, l0 = var("g, m, c, l0")

y  = dynamicsymbols('y')
yp = dynamicsymbols('y', 1)

# Kinetic energy:
T = m/2*yp*yp

# Spring energy:
# l, dl: deformed length, elongation
l = sqrt(l0*l0 + y*y)
dl = l - l0
U = c/2*dl**2 - m*g*y

# Lagrangian:
L = T - U

pprint("\nEquation of Motion:")
LM = LagrangesMethod(L, [y])
tmp = LM.form_lagranges_equations()
tmp = tmp[0]
tmp /= m
tmp = tmp.simplify()
pprint(tmp)

sub_list=[
    (y, 1 *meter),
    (c, 3 *newton/meter),
    (l0, 75*meter/100),
    ]

pprint("\nS / N:")
S = c*dl
tmp = S
tmp = tmp.subs(sub_list)
tmp /= newton
pprint(N(tmp,3))

pprint("\ntheta / °:")
theta = atan(y/l0)
tmp = theta
tmp = tmp.subs(sub_list)
tmp = deg(tmp)
pprint(N(tmp,3))

# Equation of Motion:
#                                       2
#       c⋅l₀⋅y(t)        c⋅y(t)        d
# - ────────────────── + ────── - g + ───(y(t))
#        _____________     m            2
#       ╱   2    2                    dt
#   m⋅╲╱  l₀  + y (t)
#
# S / N:
# 1.50
#
# theta / °:
# 53.1
