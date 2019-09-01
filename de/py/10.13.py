from sympy.physics.units import *
from sympy import *

prec = 3

newton = kg*m/s**2

a = 0.3 *m
qG = 25 *newton/m
g = 9.81 *m/s/s

pprint("\nmass / kg:")
G = a*qG
mass = G/g
tmp = mass
tmp /= kg
pprint(N(tmp,prec))

pprint("\nx-coord of centroid / m:")
(xS, yS) = (-2*a/3, a/2)
tmp = xS
tmp /= m
pprint(N(tmp,prec))

pprint("\ny-coord of centroid / m:")
tmp = yS
tmp /= m
pprint(N(tmp,prec))

pprint("\nIx'x' / (kg m²):")
Ix1 = mass * yS*yS
Ix2 = mass * a*a / 12
Ix3 = Ix1
Ix = Ix1 + Ix2 + Ix3
tmp = Ix
tmp /= kg*m**2
pprint(N(tmp,prec))

pprint("\nIy'y' / (kg m²):")
Iy1 = mass * a*a / 12 + (a/6)**2*mass
Iy2 = 0 + (a/3)**2*mass
Iy3 = Iy1
Iy = Iy1+Iy2+Iy3
tmp = Iy
tmp /= kg*m**2
pprint(N(tmp,prec))
