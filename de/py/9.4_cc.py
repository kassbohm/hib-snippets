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

l = 1.4 *m

theta = N(25*pi/180)

ct, st = cos(theta), sin(theta)

# Omega and Omega' from (1) and (2):
w  = Matrix([2, 0, 6]) / s
wp = Matrix([1.5, 12, 3]) / (s**2)

r = l*Matrix([0, ct, st])

vA = w.cross(r)
aA = wp.cross(r) + w.cross(w.cross(r))

pprint("\nvA / (m/s):")
tmp = vA / (m/s)
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\naA / (m/s²):")
tmp = aA / (m/(s**2))
tmp = iso_round(tmp,0.1)
pprint(tmp)

# vA / (m/s):
# ⎡-7.6⎤
# ⎢    ⎥
# ⎢-1.2⎥
# ⎢    ⎥
# ⎣2.5 ⎦
#
# aA / (m/s²):
# ⎡10.4 ⎤
# ⎢     ⎥
# ⎢-51.6⎥
# ⎢     ⎥
# ⎣-0.5 ⎦
