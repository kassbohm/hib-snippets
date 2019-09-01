from sympy.physics.units import *
from sympy import *

prec = 3

l = 3 *m
theta = 60 *pi/180

ct, st = cos(theta), sin(theta)

# Omega and Omega':
O  = Matrix([-2, -2, 0]) / s
Op = Matrix([-8, -5, -4]) / (s**2)

r = l*Matrix([0, ct, - st])

vA = O.cross(r)
aA = Op.cross(r) + O.cross(O.cross(r))

pprint("\nw'/ (1/s²):")
tmp = Op
tmp /= 1/(s**2)
pprint(N(tmp,prec))

pprint("\nw x ( w x r) / (m/s²):")
tmp = O.cross(O.cross(r))
tmp /= m/(s**2)
pprint(N(tmp,prec))

pprint("\nvA / (m/s):")
tmp = vA / (m/s)
pprint(N(tmp,prec))
pprint("\naA / (m/s²):")
tmp = aA / (m/(s**2))
pprint(N(tmp,prec))
