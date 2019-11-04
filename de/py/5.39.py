from sympy.physics.units import *
from sympy import *
from mpmath import radians

prec = 4

r = 3 * m/10
v = 2 * m/s
a = 3 * m/s**2

p1 = radians(50)
(c1, s1) = (cos(p1), sin(p1))

pprint("\nω / (1/s):")
w = v/(r*s1)
tmp = w
tmp /= 1/s
pprint(N(tmp,prec))

pprint("\nα / (1/s²):")
tmp  = (a/r - w*w*c1)/s1
tmp /= (1/s**2)
pprint(N(tmp,prec))
