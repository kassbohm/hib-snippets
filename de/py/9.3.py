from sympy.physics.units import *
from sympy import *

(w1, w2) = (0.6 /s, 8 /s)
beta = 45 *pi/180
(cb, sb) = (cos(beta), sin(beta))

w = Matrix([0, w2*cb, w1 + w2*sb])
wp = Matrix([-w1*w2*cb, 0, 0])

pprint("\nω / (1/s):")
tmp = w
pprint(N(tmp / (1/s),3))

pprint("\nω' / (1/s²):")
tmp = wp
pprint(N(tmp / (1/s/s),3))
