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

a = 0.3 *m
q = 25 *Newton/m
grav = 9.81 *m/s**2

pprint("\nMass of one part / kg:")
G = a*q
mass = G/grav
tmp = mass
tmp /= kg
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nx-position of centroid / m:")
(xS, yS) = (-2*a/3, a/2)
tmp = xS
tmp /= m
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\ny-position of centroid / m:")
tmp = yS
tmp /= m
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nIx'x' / (kg m²):")
Ix1 = mass * yS*yS
Ix2 = mass * a*a / 12
Ix3 = Ix1
Ix = Ix1 + Ix2 + Ix3
tmp = Ix
tmp /= kg*m**2
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nIy'y' / (kg m²):")
Iy1 = mass * a*a / 12 + (a/6)**2*mass
Iy2 = 0 + (a/3)**2*mass
Iy3 = Iy1
Iy = Iy1+Iy2+Iy3
tmp = Iy
tmp /= kg*m**2
tmp = iso_round(tmp,0.001)
pprint(tmp)
