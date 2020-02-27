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
        1,              #  1
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

half = S(1)/2

# ---

r = 40*cm
(a, w) = (4 /s**2, 2 /s)
beta = 30*pi/180

(cb, sb) = (cos(beta), sin(beta))
dCB = Matrix([-r*cb, r*sb, 0])

vC = Matrix([r*w, 0, 0])
aC = Matrix([r*a, 0, 0])

vBx, vBy = var("vBx, vBy")
aBx, aBy = var("aBx, aBy")

W = Matrix([0,0,-w])
A = Matrix([0,0,-a])
vB = Matrix([vBx, vBy, 0])
aB = Matrix([aBx, aBy, 0])

# Find vBx and vBy:
eq = Eq( vB, vC + W.cross(dCB) )
sol = solve(eq, [vBx, vBy])
vBx, vBy = sol[vBx], sol[vBy]
vB = Matrix([vBx, vBy, 0])

pprint("\nvB / (m/s):")
tmp = vB
tmp /= m/s
tmp = iso_round(tmp,0.01)
pprint(tmp)

# Find aBx and aBy:
eq = Eq( aB, aC + A.cross(dCB) + W.cross(W.cross(dCB)) )
sol = solve(eq, [aBx, aBy])
aBx, aBy = sol[aBx], sol[aBy]
aB = Matrix([aBx, aBy, 0])

pprint("\naB / (m/s²):")
tmp = aB
tmp /= m/s**2
tmp = iso_round(tmp,0.01)
pprint(tmp)

          # vB / (m/s):
# ⎡1.2 ⎤
# ⎢    ⎥
# ⎢0.69⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# aB / (m/s²):
# ⎡3.79⎤
# ⎢    ⎥
# ⎢0.59⎥
# ⎢    ⎥
# ⎣0.0 ⎦
