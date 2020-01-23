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

# 9.22
b, c, d, l = var("b, c, d, l")

vAz, aAz = var("vAz, aAz")

sub_list = [
    ( b    ,   2  *m       ),
    ( c    ,   1  *m       ),
    ( d    ,   S(3)/2 *m   ),
    ( vAz  , - 3  *m/s     ),
    ( aAz  ,   0  *m/s**2  ),
]

l = sqrt(b*b + c*c)

dAB = Matrix([0, b, -c])

vA = Matrix([0, 0, vAz])
aA = Matrix([0, 0, aAz])

# ω: 3 unknowns:
wx, wy, wz = var("wx, wy, wz")
w = Matrix([wx, wy, wz])

# vBe: 1 unknown:
vBe = var("vBe")
ca = b/sqrt(b*b + d*d)
sa = d/sqrt(b*b + d*d)
e = Matrix([-sa, ca, 0])
vB = vBe * e

# equations to solve:
eq1 = Eq(w.dot(dAB), 0)
eq2 = Eq(vB, vA + w.cross(dAB))
sol = solve([eq1, eq2], [wx,wy,wz,vBe])
wx, wy, wz, vBe = sol[wx], sol[wy], sol[wz], sol[vBe]
w = Matrix([wx, wy, wz])
vB = vBe * e

pprint("\n\nω / (1/s):")
tmp = w
tmp =  tmp.subs(sub_list)
tmp = tmp / (1/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\n\nvBe / (m/s):")
tmp = vBe
tmp =  tmp.subs(sub_list)
tmp = tmp / (m/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\n\nvB / (m/s):")
tmp = vB
tmp =  tmp.subs(sub_list)
tmp = tmp / (m/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)
