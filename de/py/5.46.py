from sympy.physics.units import *
from sympy import *
from mpmath import radians

l = 5 *m
vA = 6 *m/s
alpha = radians(30)
phi = radians(45)

sp, cp = sin(phi), cos(phi)
ca, ta = cos(alpha), tan(alpha)

pprint("\nω / (1/s):")
dCA = l * (cp*ta + sp)
w = vA / dCA
tmp = w
tmp /= (1/s)
pprint(N(tmp,3))

pprint("\nvB / (m/s):")
dCB = l*cp/ca
vB = w*dCB
tmp = vB
tmp /= m/s
pprint(N(tmp,3))
