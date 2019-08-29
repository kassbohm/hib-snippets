# %%
from sympy.physics.units import *
from sympy import *

t1 = 2 *s
s1 = 1 *m
v1 = 2 *m/s
t2 = 6 *s
b = 2 *m/s/s/s
a0 = 1 *m/s/s

v2 = v1 - a0*t2 + b*t2*t2/2

pprint("\nv2 / (m/s):")
pprint(v2 / (m/s))
