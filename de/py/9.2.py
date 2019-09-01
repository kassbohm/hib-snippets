from sympy.physics.units import *
from sympy import *

l = 12 *m
wx = S(6)/10 /s   # w2
wy = 0
wz = S(15)/100 /s # w1

beta = S(30)*pi/180

rx, ry, rz = 0, l * cos(beta), l * sin(beta)

wxp = S(4)/10 /s**2    # w2'
wyp = wz*wx            # as in 9.1
wzp = S(2)/10 /s**2    # w1'

# Angular velocity:
w = Matrix([wx,wy,wz])

# Position vector:
r = Matrix([rx,ry,rz])

# Velocity v:
pprint("\nv / (m/s):")
v = w.cross(r)
vval = v / (m/s)
pprint(N(vval,3))

# Acceleration a:
# Angular acceleration:
wp = Matrix([wxp,wyp,wzp])

a1 = wp.cross(r)
a2 = w.cross(w.cross(r))

a = a1 + a2
pprint("\na / (m/s²):")
aval = a / (m/s/s)
pprint(N(aval,3))
