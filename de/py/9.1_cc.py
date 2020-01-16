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

l = 12 *m
wx = S(6)/10 /s   # w2
wy = 0
wz = S(15)/100 /s # w1

beta = 30 *pi/180

# Angular velocity:
w = Matrix([wx,wy,wz])

# Position vector and angular acceleration:
(rx, ry, rz) = (0, l*cos(beta), l*sin(beta))
r = Matrix([rx,ry,rz])

# Angular acceleration:
pprint("\nw' / (1/s²):")
(wxp, wyp, wzp) = (0, wz*wx, S(8)/10 /s/s)
wp = Matrix([wxp,wyp,wzp])
tmp = wp
tmp /= (1/s**2)
pprint(tmp)

# Velocity v:
pprint("\nv / (m/s):")
v = w.cross(r)
vval = v / (m/s)
pprint(N(vval,3))

# Acceleration a:
a1 = wp.cross(r)
a2 = w.cross(w.cross(r))

a = a1 + a2
pprint("\na / (m/s²):")
aval = a / (m/s/s)
pprint(N(aval,3))
