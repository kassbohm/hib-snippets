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

la, lb, ha = var("la, lb, ha")
alpha = var("alpha")
rD, v = var("rD, v")

sub_list=[
    ( la,     2 *m       ),
    ( lb,     1 *m       ),
    ( ha,     3 *m/2     ),
    ( alpha, 45 *deg     ),
    ( rD,     5 *m  /100 ),
    ( v,     25 *m/s/100 ),
    ]

wD = v/rD

# to reproduce Hibbeler's (wrong) solution:
# wD *= 2

ca, sa = cos(alpha), sin(alpha)
t1 = Matrix([0, 0, ha])
t2 = la * Matrix([ca, 0,  sa])
t3 = lb * Matrix([sa, 0, -ca])
rC = t1 + t2 + t3
WD = Matrix([0,0,wD])

vC = WD.cross(rC)
aC = WD.cross(WD.cross(rC))

pprint("\nvC / (m/s):   ")
tmp = vC
tmp = tmp.subs(sub_list)
tmp /= m/s
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\naC / (m/s²):   ")
tmp = aC
tmp = tmp.subs(sub_list)
tmp /= m/s**2
tmp = iso_round(tmp,0.1)
pprint(tmp)

# vC / (m/s):
# ⎡0.0 ⎤
# ⎢    ⎥
# ⎢10.6⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# aC / (m/s):
# ⎡-53.0⎤
# ⎢     ⎥
# ⎢ 0.0 ⎥
# ⎢     ⎥
# ⎣ 0.0 ⎦
