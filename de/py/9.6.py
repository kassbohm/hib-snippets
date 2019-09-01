from sympy.physics.units import *
from sympy import *

mm = m/1000

ra = 200 *mm
rb = 50 *mm

wz  = 10 /s
wzp = 0 /s/s

w = Matrix([0, -4*wz, wz])
wp = Matrix([4*wz*wz, 0, 0])
r = Matrix([0, ra, rb])

pprint("\nv / (m/s):")
v = w.cross(r)
tmp = v
pprint(N(tmp / (m/s),3))

pprint("\na / (m/s²):")
a = wp.cross(r) + w.cross(v)
tmp = a
pprint(N(tmp / (m/s/s),3))

pprint("\n|ω| / (1/s):")
tmp = w.norm()
pprint(N(tmp / (1/s),3))

pprint("\n|ω'| / (1/s²):")
tmp = wp.norm()
pprint(N(tmp / (1/s/s),3))
