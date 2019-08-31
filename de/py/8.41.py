from sympy.physics.units import *
from sympy import *

(mAB, mS) = (2 *kg, 4 *kg)
(r, l) = (S(15)/100 *m, S(15)/10 *m)
w1 = 5 *1/s

LC1 = mS *r*r *w1

J2 = mAB/12 * l*l + 2*mS*(l/2)**2

pprint("\nw2 / (1/s):")
w2 = LC1/J2
tmp = w2
tmp /= 1/s
pprint(N(tmp,3))
