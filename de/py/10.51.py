from sympy.physics.units import *
from sympy import *

prec = 3
g  = S(981)/100 *m/s**2

newton = kg*m/s**2
Nm = newton*m

b  =  1 *m
M  = 50 *Nm
w  = 10 /s
mu = 5 *kg/m

wx, wy, wz = var("omega_x, omega_y, omega_z")
ax, ay, az = var("alpha_x, alpha_y, alpha_z")
Mx, My, Mz = var("M_x, M_y, M_z")
Ixx, Iyy, Izz = var("I_xx, I_yy, I_zz")
Ixy, Ixz, Iyz = var("I_xy, I_xz, I_yz")

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

# Substitutions:
# --- 1 ---
sub_list = [ (wx, 0), (wy, 0), (wz, w), (ax, 0), (ay, 0)]
eq1 = eq1.subs(sub_list)
eq2 = eq2.subs(sub_list)
eq3 = eq3.subs(sub_list)

# --- 2 ---
Bx, By = var("Bx, By")
M_x = - 3*b*By - mu*b*g*3*b/2
M_y = 3*b*Bx
M_z = M
eq1 = eq1.subs(Mx, M_x)
eq2 = eq2.subs(My, M_y)
eq3 = eq3.subs(Mz, M_z)

# --- 3 ---
I_xy = 0
I_xz = 0
I_yz = - b*b/2*mu*b - b*3*b/2*mu*b
I_zz = b*b/3 * mu*b + b**2*mu*b

sub_list = [ (Ixy, I_xy), (Ixz, I_xz), (Iyz, I_yz), (Izz, I_zz)]
eq1 = eq1.subs(sub_list)
eq2 = eq2.subs(sub_list)
eq3 = eq3.subs(sub_list)

sol = solve([eq1, eq2, eq3], [Bx, By, az])
Bx = sol[Bx]
By = sol[By]
az = sol[az]
pprint("\nBx / N:")
tmp = Bx
tmp /= newton
pprint(N(tmp,prec))

pprint("\nBy / N:")
tmp = By
tmp /= newton
pprint(N(tmp,prec))

pprint("\nw'z / (1/s²):")
tmp = az
tmp /= 1/s**2
pprint(N(tmp,prec))

# Bx / N:
# -25.0
#
# By / N:
# -358.
#
# w'z / (1/s²):
# 7.50
