from sympy import *
from sympy.physics.units import newton, meter, second, kilogram

newton = kg*m/s**2
g = 981 *m/s**2 / 100

G = 8 * newton
r = S(3)/10 *meter
l = S(5)/10 *meter
ws = 300 / second

theta = S(40)*pi/180
# Solve numerically:
theta = N(theta)

wp = var("omega_p")

c, s = cos(theta), sin(theta)

m = G/g

J = m*r*r/4 + l*l*m
Jz = m*r*r/2

tmp = J / kilogram/meter/meter
pprint("\nJ / (kg m²): ")
pprint(N(tmp,3))

# Terms in equation divided by sin(theta):
t1 = G * l
t2 = - J *wp*wp *c
t3 = Jz * wp *(wp*c + ws)

sol = solve(t1 - (t2 + t3))
pprint("\nSolutions for wp:")
pprint(sol)

# J / (kg m²):
# 0.222
#
# Solutions for wp:
# ⎡0.365053676362308  77.098675580718⎤
# ⎢─────────────────, ───────────────⎥
# ⎣      second            second    ⎦
