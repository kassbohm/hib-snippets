from sympy.physics.units import *
from sympy import *

prec = 3

g = S(981)/100 *m/s**2

newton = kg*m/s**2
kN = 1000 * newton
Nm = newton *m
kNm = 1000 *Nm

G = 50 *kN
k = S(2)/10 *m
wy = 2 /s
wS = 100 /s

pprint("\nm / kg:")
mass = G / g
tmp = mass
tmp /= kg
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nIzz / (kg m²):")
mass = G / g
Izz = mass*k**2
tmp = Izz
tmp /= kg *m**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nMx / (kNm):")
Mx = Izz*wy*wS
tmp = Mx
tmp /= kNm
tmp = N(tmp,prec)
pprint(tmp)
pprint("\nMx / (Nm):")
tmp = Mx
tmp /= Nm
tmp = N(tmp,prec)
pprint(tmp)

# m / kg:
# 5.10e+3
#
# Izz / (kg m²):
# 204.
#
# Mx / (kNm):
# 40.8
#
# Mx / (Nm):
# 4.08e+4
