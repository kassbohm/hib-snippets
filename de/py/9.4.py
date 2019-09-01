from sympy.physics.units import *
from sympy import *

prec = 3

l = 1.4 *m
theta = N(25*pi/180)

ct, st = cos(theta), sin(theta)

# Omega and Omega' from (1) and (2):
w  = Matrix([2, 0, 6]) / s
wp = Matrix([1.5, 12, 3]) / (s**2)


r = l*Matrix([0, ct, st])

vA = w.cross(r)
aA = wp.cross(r) + w.cross(w.cross(r))

pprint("\nvA / (m/s):")
tmp = vA / (m/s)
pprint(N(tmp,prec))
pprint("\naA / (m/s²):")
tmp = aA / (m/(s**2))
pprint(N(tmp,prec))
