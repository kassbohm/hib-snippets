from sympy.physics.units import *
from sympy import *

# solution as in book - except definition of alpha:
# here: alpha between z and rod
# Hibbeler: alpha between y and rod

# Given:
# r can also be computed from l and alpha, but:
l, alpha, m, r = var("l, alpha, m, r")

# phi and phi':
phi, phi1 = var("phi, phi1")

ca, sa = cos(alpha), sin(alpha)

# kinetic energy:
# translation:
vMx = l*phi1*sa*ca
T1 = m*vMx*vMx/2

# rotation:
w = phi1*sa
wi = 0
wj = w*sa
wk = - w*ca

J = m*r*r
Ji = J/4
Jj = J/2
Jk = J/4

T2 = ( Ji*wi*wi + Jj*wj*wj + Jk*wk*wk ) / 2

# total:
T = T1 + T2
pprint(T)

dTdp1 = diff(T,phi1)
tmp = dTdp1/phi1/m
tmp = tmp.simplify()
pprint(tmp)
