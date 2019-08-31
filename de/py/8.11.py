from sympy.physics.units import *
from sympy import *

mm = m/1000
newton=kg*m/s**2
kN = 1000*newton

g = S(981)/100 *m/s/s
G = 85 *kN

r = 375 *mm
k = S(14)/10 *m
t1 = 5 *s
v0 = 360 *m/s
F1 = 25 *kN
F2 = 4 *kN

mass = G/g
J = mass*k*k

v1 = v0 + (F1 + F2)/mass * t1
w1 = r *(F1 - F2)/J *t1
pprint("\nv / (m/s):")
pprint(N(v1/ (m/s),3))
pprint("\nω / (1/s):")
pprint(N(w1 / (1/s),3))
