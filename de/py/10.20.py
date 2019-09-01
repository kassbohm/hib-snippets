from sympy.physics.units import *
from sympy import *

prec = 3

Ixx, Iyy, Izz, Ixy, Ixz, Iyz = var("Ixx, Iyy, Izz, Ixy, Ixz, Iyz")

ux, uy, uz = var("ux, uy, uz")

I = Matrix([
    [Ixx, Ixy, Ixz],
    [Ixy, Iyy, Iyz],
    [Ixz, Iyz, Izz],
    ])

I_Hib = Matrix([
    [ Ixx, -Ixy, -Ixz],
    [-Ixy,  Iyy, -Iyz],
    [-Ixz, -Iyz,  Izz],
    ])

u = Matrix([ux, uy, uz])
uT = u.transpose()

Iuu = uT.dot(I.dot(u))
pprint(Iuu.expand())

#       2                                     2                       2
# Ixx⋅ux  + 2⋅Ixy⋅ux⋅uy + 2⋅Ixz⋅ux⋅uz + Iyy⋅uy  + 2⋅Iyz⋅uy⋅uz + Izz⋅uz

# Hib:
Iuu_Hib = uT.dot(I_Hib.dot(u))
pprint(Iuu_Hib.expand())

#       2                                     2                       2
# Ixx⋅ux  - 2⋅Ixy⋅ux⋅uy - 2⋅Ixz⋅ux⋅uz + Iyy⋅uy  - 2⋅Iyz⋅uy⋅uz + Izz⋅uz

# ---

pprint("\n")

a = S(6)/10 *m
h = S(4)/10 *m
mb = 4 *kg/m

m1 = 2*a*mb
m2 =   a*mb
m3 =   h*mb

Ixx1 = m1/12*2*a*2*a
Ixx1 += a*a*m1
Ixx2 = m2 * 2*a*2*a
Ixx3 = m3/12*h*h
Ixx3 += 2*a*2*a*m3
Ixx3 += h*h*m3/4
I_xx = Ixx1 + Ixx2 + Ixx3

Iyy1 = 0
Iyy2 = m2/12*a*a
Iyy2 += m2*a*a/4
Iyy3 = m3/12*h*h
Iyy3 += m3*a*a
Iyy3 += m3*h*h/4
I_yy = Iyy1 + Iyy2 + Iyy3

Izz1 = Ixx1
Izz2 = m2/12 *a*a
Izz2 += m2*a*a/4
Izz2 += m2*2*a*2*a
Izz3 = m3*(a*a + 2*a*2*a)
I_zz = Izz1 + Izz2 + Izz3

Ixy1 = 0
Ixy2 = m2*(- a/2 * 2*a)
Ixy3 = m3*(- a * 2*a)
I_xy = Ixy1 + Ixy2 + Ixy3

Iyz1 = 0
Iyz2 = 0
Iyz3 = m3*(- 2*a * h/2)
I_yz = Iyz1 + Iyz2 + Iyz3

Ixz1 = 0
Ixz2 = 0
Ixz3 = m3*(- a * h/2)
I_xz = Ixz1 + Ixz2 + Ixz3

d0a = Matrix([a, 2*a, h])
u = d0a / d0a.norm()
u_x, u_y, u_z = u[0], u[1], u[2]

sub_list=[
    (Ixx, I_xx),
    (Iyy, I_yy),
    (Izz, I_zz),
    (Ixy, I_xy),
    (Iyz, I_yz),
    (Ixz, I_xz),
    (ux, u_x),
    (uy, u_y),
    (uz, u_z),
    ]

pprint("\nIuu / (kg m²):")
tmp = Iuu.subs(sub_list)
tmp /= kg * m**2
pprint(N(tmp,prec))
