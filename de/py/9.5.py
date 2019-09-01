from sympy.physics.units import *
from sympy import *

# Solution in book incorrect.

mm = m / 1000

ra = 160 *mm
rb = 80 *mm

wz  = 5 /s
wzp = 2 /s/s

w = Matrix([0, -2*wz, wz])
wp = Matrix([2*wz*wz, -2*wzp, wzp])
r = Matrix([0, ra, rb])

pprint("\nv / (m/s):")
v = w.cross(r)
tmp = v
pprint(N(tmp / (m/s),3))

pprint("\na / (m/s²):")
a = wp.cross(r) + w.cross(v)
tmp = a
pprint(N(tmp / (m/s/s),3))
