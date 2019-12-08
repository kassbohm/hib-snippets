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

prec = 3

(r, rA, l) = (3*cm, 2*cm, 8*cm)
(w, a) = (6 /s, 12 /s**2)
beta = 60*pi/180
(cb, sb) = (cos(beta), sin(beta))

w1 = Matrix([0,0,w])
a1 = Matrix([0,0,a])

vAx, vAy, aAx, aAy = var("vAx, vAy, aAx, aAy")
vBx, aBx = var("vBx, aBx")

vA = Matrix([vAx, vAy, 0])
aA = Matrix([aAx, aAy, 0])

vB = Matrix([vBx, 0, 0])
aB = Matrix([aBx, 0, 0])

aP = Matrix([0, r*w**2, 0])

h = r - rA
dPA = Matrix([0, h, 0])
dAB = Matrix([l*cb, l*sb, 0])

# Rigid body 1:
# Find vAx and vAy:
eq = Eq( vA, w1.cross(dPA) )
sol = solve(eq, [vAx, vAy])
vAx, vAy = sol[vAx], sol[vAy]
vA = Matrix([vAx, vAy, 0])

pprint("\nvA / (cm/s):")
tmp = vA
tmp /= cm/s
tmp = N(tmp,prec)
pprint(tmp)

# Find aAx and aAy:
eq = Eq( aA, aP + a1.cross(dPA) + w1.cross(vA) )
sol = solve(eq, [aAx, aAy])
aAx, aAy = sol[aAx], sol[aAy]
aA = Matrix([aAx, aAy, 0])

pprint("\naA / (cm/s²):")
tmp = aA
tmp /= cm/s**2
tmp = N(tmp,prec)
pprint(tmp)

# Rigid body 2:
w2, a2 = var("w2, a2")
W2 = Matrix([0,0,w2])
A2 = Matrix([0,0,a2])

# Find vBx and w2:
eq = Eq( vB, vA + W2.cross(dAB) )
sol = solve(eq, [vBx, w2])
vBx, w2 = sol[vBx], sol[w2]
vB = Matrix([vBx, 0, 0])

pprint("\nvB / (cm/s):")
tmp = vB
tmp /= cm/s
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nw2 / (1/s):")
W2 = Matrix([0,0,w2])
tmp = W2
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)

# Find aBx and a2:
eq = Eq( aB, aA + A2.cross(dAB) + W2.cross(W2.cross(dAB)) )
sol = solve(eq, [aBx, a2])
aBx, a2 = sol[aBx], sol[a2]

pprint("\naB / (cm/s²):")
aB = Matrix([aBx, 0, 0])
tmp = aB
tmp /= cm/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\na2 / (1/s²):")
A2 = Matrix([0,0,a2])
tmp = A2
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)

# vA / (cm/s):
# ⎡-6.0⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# aA / (cm/s²):
# ⎡-12.0⎤
# ⎢     ⎥
# ⎢72.0 ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# vB / (cm/s):
# ⎡-6.0⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# w2 / (1/s):
# ⎡0⎤
# ⎢ ⎥
# ⎢0⎥
# ⎢ ⎥
# ⎣0⎦
#
# aB / (cm/s²):
# ⎡113.0⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# a2 / (1/s²):
# ⎡  0  ⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣-18.0⎦
