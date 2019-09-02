from sympy.physics.units import *
from sympy import *

prec = 3
g  = S(981)/100 *m/s**2
newton = kg*m/s**2
Nm = newton*m

wx, wy, wz = var("omega_x, omega_y, omega_z")
ax, ay, az = var("alpha_x, alpha_y, alpha_z")
Mx, My, Mz = var("Mx, My, Mz")
Ixx, Iyy, Izz = var("Ixx, Iyy, Izz")
Ixy, Ixz, Iyz = var("Ixy, Ixz, Iyz")

# Euler:
eq1 = Eq( Mx, Ixx*ax
    - (Iyy - Izz) * wy*wz
    + Ixy * (ay - wz*wx)
    + Iyz * (wy**2 - wz**2)
    + Ixz * (az + wx*wy)
    )
eq2 = Eq( My, Iyy*ay
    - (Izz - Ixx) * wz*wx
    + Iyz * (az - wx*wy)
    + Ixz * (wz**2 - wx**2)
    + Ixy * (ax + wy*wz)
    )
eq3 = Eq( Mz, Izz*az
    - (Ixx - Iyy) * wx*wy
    + Ixz * (ax - wy*wz)
    + Ixy * (wx**2 - wy**2)
    + Iyz * (ay + wz*wx)
    )

mass, b, d, w, a = var("m, b, d, omega, alpha")

Ixx = mass*b**2/12
Iyy = 0
Izz = mass*b**2/12
(ax, ay, az) = (0, (d + b/2) * w**2, b/2 * a)
(wx, wy, wz) = (0, 0, w)
(alphax, alphay, alphaz) = (-a, 0, 0)

sub_list = [
    (mass, S(25)/10 *kg),
    (b, S(9)/10 *m),
    (d, S(6)/10 *m),
    (w, 3 /s),
    (a, 2 /s**2),
    ]

pprint("\nIxx, Izz / (kg m²):")
tmp = Ixx.subs(sub_list)
tmp /= kg*m**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nay / (m/s²):")
tmp = ay.subs(sub_list)
tmp /= m/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\naz / (m/s²):")
tmp = az.subs(sub_list)
tmp /= m/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nalphax / (1/s²):")
tmp = alphax.subs(sub_list)
tmp /= 1/s**2
tmp = N(tmp,prec)
pprint(tmp)

# Ixx, Izz / (kg m²):
# 0.169
#
# ay / (m/s²):
# 9.45
#
# az / (m/s²):
# 0.900
#
# alphax / (1/s²):
# -2.00

eq1 = eq1.subs(sub_list)
eq2 = eq2.subs(sub_list)
eq3 = eq3.subs(sub_list)
