from sympy.physics.units import *
from sympy import *

prec = 4
newton = kg * m / s**2
Nm = newton*m
mm = m / 1000

# ---

mass = S(8)/10 *kg
l = 150 *mm
w = 2 / s

theta = var("theta")
C, S = cos(theta), sin(theta)

# theta':
tp = 6 / s**2

wx = w*S
wy = w*C
wz = tp

ax = 2*tp*C
ay = -2*tp*S
az = 0

Mz = var("M_z")
Ixx, Iyy, Izz = var("Ixx, Iyy, Izz")
Ixy, Ixz, Iyz = var("Ixy, Ixz, Iyz")

# Euler:
eq3 = Eq( Mz, Izz*az
    - (Ixx - Iyy) * wx*wy
    + Ixz * (ax - wy*wz)
    + Ixy * (wx**2 - wy**2)
    + Iyz * (ay + wz*wx)
    )

eq3 = eq3.subs([
    (Ixy, 0),
    (Ixz, 0),
    (Iyz, 0),
    (Ixx,0),
    (Iyy,mass*l**2/12),
    (Izz,mass*l**2/12),
    ])

pprint(eq3)

pprint("\nMz / (Nm):")
Mz = solve(eq3, Mz)[0]
tmp = Mz / Nm
pprint(tmp)

#                       2
#       3⋅kilogram⋅meter ⋅sin(θ)⋅cos(θ)
# M_z = ───────────────────────────────
#                           2
#                 500⋅second
#
# Mz / (Nm):
# 3⋅sin(2⋅θ)
# ──────────
#    1000
