from sympy.physics.units import *
from sympy import *

l = 2 *m
beta = N(20 *pi/180)
wz = 4 /s
wzp = 3 /s/s

cb, sb = cos(beta), sin(beta)
c2b, s2b = cos(2*beta), sin(2*beta)

w = Matrix([0, -wz*cb/sb, 0])

pprint("\nw / (1/s):")
tmp = w
tmp /= 1/s
pprint(N(tmp,3))

rA = l * Matrix([0, c2b, s2b])
pprint("\nrA / m:")
tmp = rA / m
pprint(N(tmp,3))

vA = w.cross(rA)

pprint("\nvA / (m/s):")
tmp = vA / (m/s)
pprint(N(tmp,3))

wp_xyz = Matrix([0, -wzp*cb/sb, 0])
wp_e = Matrix([wz*wz*cb/sb, 0, 0])

wp = wp_xyz + wp_e
pprint("\nwp / (1/s²):")
tmp = wp / (1/s/s)
pprint(N(tmp,3))

aA = wp.cross(rA) + w.cross(vA)
pprint("\naA / (m/s²):")
tmp = aA / (m/s/s)
pprint(N(tmp,3))
