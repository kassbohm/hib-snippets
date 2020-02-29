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

grav = 981 *m/s**2 / 100

mass, a, mu = var("mass, a, mu")

sub_list = [
    (mass, 800 *kg),
    (a, S(5)/10 *m/s**2),
    (mu, S(1)/10),
    ]

G = mass*grav
c = sqrt(2)/2

Nv, S = var("Nv, S")
eq1 = Eq(mass*a, S*c - mu*Nv)
eq2 = Eq(G, Nv + S*c)
sol = solve([eq1, eq2],[Nv, S])

pprint("\nN / Newton:")
Nv = sol[Nv]
tmp = Nv
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= Newton
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\nS / Newton:")
S = sol[S]
tmp = S
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= Newton
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\nphi / deg:")
phi = asin(mu * Nv / S)
tmp = phi
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp = tmp * 180 / pi
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\ntheta / deg:")
tmp = 45 - tmp
tmp = iso_round(tmp,0.001)
pprint(tmp)

# N / Newton:
# 6770.9
#
# S / Newton:
# 1523.2
#
# phi / deg:
# 26.392
#
# theta / deg:
# 18.608
