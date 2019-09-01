from sympy.physics.units import s, m, kg
from sympy import *

# 9.21 = 10.24 / 10.25

b, c, d, l = var("b, c, d, l")

cm = m/100
newton = kg*m/s**2

sub_list = [
    (b, 2 *cm),
    (c, 3 *cm),
    (d, 6 *cm),
    (l, 7 *cm),
]

vAz = - 8 *cm/s
aAz = - 5 *cm/s/s
vA = Matrix([0, 0, vAz])
aA = Matrix([0, 0, aAz])

# Eight unknowns:
# ω and α:
wx, wy, wz = var("wx, wy, wz")
ax, ay, az = var("ax, ay, az")
vBx, aBx = var("vBx, aBx")

# Vectors in which unknowns appear:
w = Matrix([wx, wy, wz])
a = Matrix([ax, ay, az])
vB = Matrix([vBx, 0, 0])
aB = Matrix([aBx, 0, 0])

dAB = Matrix([b, d, -c])

# equations to solve:
eq1 = Eq(w.dot(dAB), 0)
eq2 = Eq(vB, vA + w.cross(dAB))
sol = solve([eq1, eq2], [wx,wy,wz,vBx])

wx, wy, wz, vBx = sol[wx], sol[wy], sol[wz], sol[vBx]
w = Matrix([wx, wy, wz])

vB = Matrix([vBx, 0, 0])

pprint("\n\nω / (1/s):")
tmp =  w.subs(sub_list)
tmp = tmp / (1/s)
pprint(N(tmp,3))

pprint("\n\nvB / (cm/s):")
tmp = vB.subs(sub_list)
tmp = tmp / (cm/s)
pprint(N(tmp,3))

# equations to solve:
eq1 = Eq(aB, aA + a.cross(dAB) + w.cross(w.cross(dAB)))
# eq2 = Eq(a.dot(dAB) + w.dot(w.cross(dAB)), 0)
# Note, that: w . ( w x d ) = 0, so that:
eq2 = Eq(a.dot(dAB), 0)

sol = solve([eq1, eq2], [ax,ay,az,aBx])
ax, ay, az, aBx = sol[ax], sol[az], sol[az], sol[aBx]
a = Matrix([ax, ay, az])
aB = Matrix([aBx, 0, 0])

pprint("\n\naB / (cm/s²):")
tmp = aB.subs(sub_list)
tmp = tmp / (cm/s/s)
pprint(N(tmp,3))

# 10.24/10.25:
# xB, yB, zA = b, d, c

pprint("\n\nvG / (cm/s):")
vG = vA + w.cross(dAB)/2
tmp = vG.subs(sub_list)
tmp /= (cm/s)
tmp = N(tmp,3)
pprint(tmp)

# ω²:
v2 = vG.dot(vG)
w2 = w.dot(w)

# mass and moment of inertia (component z'z'):
mass = 6 *kg
Izz = mass * l**2/12

pprint("\n\nT / (Nm):")
T = (mass*v2 + Izz*w2)/2
tmp = T.subs(sub_list)
tmp /= (newton*m)
tmp = N(tmp,3)
pprint(tmp)

pprint("\n\nTe / (Nm):")
g = 9.81 *m/s**2
dW = 0.5*mass*g*c
Te = T + dW
tmp = Te.subs(sub_list)
tmp /= (newton*m)
tmp = N(tmp,3)
pprint(tmp)
