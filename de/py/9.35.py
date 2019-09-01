from sympy.physics.units import *
from sympy import *

prec = 3

l, c = 12 *m, S(15)/10 *m
t1, w1 = 60 *pi/180, S(5)/10 / second
wz, wzp = S(75)/100 / second, 2 / second**2

c1, s1 = cos(t1), sin(t1)

# Angular velocity of (x,y,z)-system:
Oz = Matrix([0, 0, wz])
Ozp= Matrix([0, 0, wzp])

# Angular velocity of AB relative to (x,y,z):
Oy = Matrix([0, -w1, 0])

# Position of A:
rA = Matrix([c, 0, 0])

# Vector from A to B:
dAB = Matrix([l*c1, 0, l*s1])

# Position of B:
rB = rA + dAB

# ----------------
# --- Velocity ---
# ----------------

# Velocity of A relative to inertial frame:
vA = Oz.cross(rA)
tmp = vA
tmp /= m/second
pprint("\nvA / (m/s):")
pprint(N(tmp, prec))

# Velocity of B relative to (x,y,z):
vBxyz = Oy.cross(dAB)
tmp = vBxyz
tmp /= m/second
pprint("\nvBxyz / (m/s):")
pprint(N(tmp, prec))

# Velocity of B relative to inertial frame:
vB = vA + Oz.cross(dAB) + vBxyz
tmp = vB
tmp /= m/second
pprint("\nvB / (m/s):")
pprint(N(tmp, prec))

# --------------------
# --- Acceleration ---
# --------------------

# Acceleration of A relative to inertial frame:
aA = Ozp.cross(rA) + Oz.cross(Oz.cross(rA))
tmp = aA
tmp /= m/second**2
pprint("\naA / (m/s²):")
pprint(N(tmp, prec))

# Acceleration of B relative to (x,y,z):
aBxyz = Oy.cross(Oy.cross(dAB))
tmp = aBxyz
tmp /= m/second**2
pprint("\naBxyz / (m/s²):")
pprint(N(tmp, prec))

# Acceleration of B relative to inertial frame:
aB = aA + Ozp.cross(dAB) + Oz.cross(Oz.cross(dAB)) + 2*Oz.cross(vBxyz) + aBxyz
tmp = aB
tmp /= m/second**2
pprint("\naB / (m/s²):")
pprint(N(tmp, prec))
