from sympy.physics.units import *
from sympy import *

# Given quantities:
l = 12 *m
wx = S(6)/10 /s   # w2
wy = 0
wz = S(15)/100 /s # w1

beta = 30 *pi/180

# Angular velocity:
w = Matrix([wx,wy,wz])

# Position vector and angular acceleration:
(rx, ry, rz) = (0, l*cos(beta), l*sin(beta))
(wxp, wyp, wzp) = (0, wz*wx, S(8)/10 /s/s)
r = Matrix([rx,ry,rz])
wp = Matrix([wxp,wyp,wzp])

# Velocity v:
pprint("\nv / (m/s):")
v = w.cross(r)
vval = v / (m/s)
pprint(N(vval,3))

# Acceleration a:
a1 = wp.cross(r)
a2 = w.cross(w.cross(r))

a = a1 + a2
pprint("\na / (m/s²):")
aval = a / (m/s/s)
pprint(N(aval,3))
