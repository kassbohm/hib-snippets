from sympy.physics.units import *
from sympy import *

# x-component of acceleration is wrong in Hibbeler's book

d = 1 *m
theta = 30 *pi/180

(w1, w2)   = (3 /s, S(3)/2 /s)
(w1p, w2p) = (2 /s**2, 4 /s**2)

# ---
(ct, st) = (cos(theta), sin(theta))

w = Matrix([w2, 0, w1])
wp = Matrix([w2p, w1*w2, w1p])
r = Matrix([0, d*ct, d*st])

pprint("\nv / (m/s):")
v = w.cross(r)
tmp = v
tmp /= m/s
pprint(N(tmp,3))

pprint("\na / (m/s²):")
a = wp.cross(r) + w.cross(v)
tmp = a
tmp /= m/s**2
pprint(N(tmp,3))
