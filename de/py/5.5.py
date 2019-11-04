from sympy.physics.units import *
from sympy import *

prec = 3

w0 = 50 / s
b = S(6)/100 /s**2

dt = 10*2*pi
w2 = w0*w0 + 2*b*dt**3/3

pprint("\nw / (1/s):")
tmp = sqrt(w2)
tmp /= 1/s
tmp = N(tmp, prec)
pprint(tmp)

# w / (1/s):
# 111.
