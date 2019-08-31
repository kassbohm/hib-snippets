# %%
from sympy.physics.units import *
from sympy import *

t1 = 2 *s
s1 = 1 *m
v1 = 2 *m/s
t2 = 6 *s
b  = 2 *m/s/s/s
a0 = 1 *m/s/s

s2 = s1 + v1*t2 - a0*t2*t2/2 +b*t2*t2*t2/6
v2 = v1 - a0*t2 + b*t2*t2/2

pprint("\ns2 / m:")
tmp = s2
tmp /= m
pprint(N(tmp,3))

pprint("\nv2 / (m/s):")
tmp = v2
tmp /= m/s
pprint(N(tmp,3))
