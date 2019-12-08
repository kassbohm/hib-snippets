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

rB, l = 150 *mm, 500 *mm
w, a  = 8 /s, 16 /s**2
phi   = 30 *deg
w = Matrix([0,0,-w])
a = Matrix([0,0,-a])


cp, sp = cos(phi), sin(phi)
er = Matrix([ -cp, sp, 0 ])
el = Matrix([  sp, cp, 0 ])
dBA = l*el

vAx = var("vAx")
vA = Matrix([vAx, 0, 0])

vBx, vBy = var("vBx, vBy")
vB = Matrix([vBx, vBy, 0])

wS = var("wS")
WS = Matrix([0,0,-wS])

# ---

eq = Eq(vB, w.cross(rB*er))
sol = solve(eq, [vBx, vBy])
vBx, vBy = sol[vBx], sol[vBy]
vB = Matrix([vBx, vBy, 0])

pprint("\nvB / (m/s):")
tmp = vB
tmp /= (m/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)

eq = Eq(vA, vB + WS.cross(dBA))
sol = solve(eq, [vAx, wS])
vAx, wS = sol[vAx], sol[wS]
vA = Matrix([vAx, 0, 0])
WS = Matrix([0,0,-wS])

pprint("\nvA / (m/s):")
tmp = vA
tmp /= (m/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)

pprint("\nWS / (1/s):")
tmp = WS
tmp /= (1/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)
           
# vB / (m/s):
# ⎡0.6 ⎤
# ⎢    ⎥
# ⎢1.04⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# vA / (m/s):
# ⎡2.4⎤
# ⎢   ⎥
# ⎢0.0⎥
# ⎢   ⎥
# ⎣0.0⎦
#
# WS / (1/s):
# ⎡ 0.0 ⎤
# ⎢     ⎥
# ⎢ 0.0 ⎥
# ⎢     ⎥
# ⎣-4.16⎦
