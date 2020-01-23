from sympy.physics.units import *
from sympy import *

# Rounding:
import decimal
from decimal import Decimal as DX
def iso_round(obj, pv, rounding=decimal.ROUND_HALF_EVEN):
    import sympy
    """
    Rounding acc. to DIN EN ISO 80000-1:2013-08
    place value = Rundestellenwert
    """
    assert pv in set([
        # place value   #  round to:
        100,            #  3rd last digit before decimal
        10,             #  2nd last
        1,              #  last
        0.1,            #  1st digit after decimal
        0.01,           #  2nd
        0.001,          #  3rd
        0.0001,         #  4th
        0.00001,        #  5th
        0.000001,       #  6th
        0.0000001,      #  7th
        0.00000001,     #  8th
        0.000000001,    #  9th
        0.0000000001,   # 10th
        ])
    try:
        tmp = DX(str(float(obj)))
        obj = tmp.quantize(DX(str(pv)), rounding=rounding)
    except:
        for i in range(len(obj)):
            tmp = DX(str(float(obj[i])))
            obj[i] = tmp.quantize(DX(str(pv)), rounding=rounding)
    return obj

# LateX:
kwargs = {}
kwargs["mat_str"] = "bmatrix"
kwargs["mat_delim"] = ""
# kwargs["symbol_names"] = {FB: "F^{\mathsf B}", }

# Units:
(k, M, G ) = ( 10**3, 10**6, 10**9 )
(mm, cm, deg) = ( m/1000, m/100, pi/180)
Newton = kg*m/s**2
Pa     = Newton/m**2
MPa    = M*Pa
GPa    = G*Pa
kN     = k*Newton

# ---

#  9.21, 10.24 and 10.25
b, c, d, l = var("b, c, d, l")

vAz, aAz = var("vAz, aAz")

sub_list = [
    ( b    ,   2  *cm      ),
    ( c    ,   3  *cm      ),
    ( d    ,   6  *cm      ),
    ( vAz  , - 8  *cm/s    ),
    ( aAz  , - 5  *cm/s**2 ),
]

l = sqrt(b*b + c*c + d*d)

# ---

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
tmp = w
tmp =  tmp.subs(sub_list)
tmp = tmp / (1/s)
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\n\nvB / (cm/s):")
tmp = vB
tmp =  tmp.subs(sub_list)
tmp = tmp / (cm/s)
tmp = iso_round(tmp,0.1)
pprint(tmp)

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
tmp = aB
tmp =  tmp.subs(sub_list)
tmp = tmp / (cm/s/s)
tmp = iso_round(tmp,0.1)
pprint(tmp)

# 10.24/10.25:
# xB, yB, zA = b, d, c

pprint("\n\nvG / (cm/s):")
vG = vA + w.cross(dAB)/2
tmp = vG.subs(sub_list)
tmp /= (cm/s)
tmp = iso_round(tmp,0.1)
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
tmp /= (Newton*m)
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\n\nTe / (Nm):")
g = 9.81 *m/s**2
dW = 0.5*mass*g*c
Te = T + dW
tmp = Te.subs(sub_list)
tmp /= (Newton*m)
tmp = iso_round(tmp,0.001)
pprint(tmp)
