from sympy.physics.units import *
from sympy import *

g = 981 *m/s**2 /100
mm = m/1000

F = 40 *newton
G = 1250 *newton
k = 240 *mm
t1 = 3 *s
r = 375 *mm

mass = G/g
J = mass*k*k

w1 = r * F/J *t1
pprint("\nω / (1/s):")
pprint(N(w1 / (1/s),3))
