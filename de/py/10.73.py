from sympy.physics.units import *
from sympy import *

prec = 3
g  = S(981)/100 *m/s**2
grams = kg/1000
mm = m/1000

mass = 200 *grams
d = 30 *mm
k = 20 *mm
wP = 2 /s

Jzz = mass*k**2

wR = mass*g*d/(Jzz*wP)

pprint("\nwR / (1/s):")
tmp = wR
tmp /= (1/s)
tmp = N(tmp,prec)
pprint(tmp)

# wR / (1/s):
# 368.
