from sympy.physics.units import *
from sympy import *

# precision:
prec = 3

r = 6 *m
theta = 45 *pi/180
w1, w2 = 6 / second, 4 / second
v, a = 5 * m/second, 8 *m/second**2

c, s = cos(theta), sin(theta)

# Angular velocity of (x,y,z)-system:
Oz = Matrix([0, 0, w1])  # constant direction
Ox = Matrix([w2, 0, 0])  # direction changing with time
O = Oz + Ox   # direction changing with time

# Position of P:
dOP = Matrix([0, r*c, r*s])

# ----------------
# --- Velocity ---
# ----------------

# Velocity of P relative to (x,y,z)-system:
vPxyz = Matrix([0, v*c, v*s])

vP = O.cross(dOP) + vPxyz
tmp = vP
tmp /= m/second
pprint("\nvP / (m/s):")
pprint(N(tmp, prec))

# --------------------
# --- Acceleration ---
# --------------------

# Angular acceleration of (x,y.z)-system due to the
# change of *direction* of vector of angular velocity:
Op = Oz.cross(Ox)

# Acceleration of P relative to (x,y,z)-system:
aPxyz = Matrix([0, a*c, a*s])

aP = Op.cross(dOP) + O.cross(O.cross(dOP)) + 2 * O.cross(vPxyz) + aPxyz
tmp = aP
tmp /= m/second**2
pprint("\naP / (m/s²):")
pprint(N(tmp, prec))
