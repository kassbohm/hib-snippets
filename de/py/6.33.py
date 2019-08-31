from sympy.physics.units import *
from sympy import *

prec = 5

newton = kg*m/s**2
kN = 1000*newton

g = S(981)/100 *m/s**2

qG, d, b, r, a = var("qG, d, b, r, a")

sub_list = [
    (qG, S(75)/10 *kN/m),
    (d, 2 *m),
    (b, 5 *m),
    (r, 1 *m),
    (a, S(1)/2 *m/s**2),
    ]

mu = qG/g

# cos(alpha):
ca = sqrt(2)/2

pprint("\nR / N:")
R = mu*(g+a)*2*b
tmp = R
tmp = tmp.subs(sub_list)
tmp /= newton
pprint(N(tmp,prec))

pprint("\nM / Nm:")
M = R * r
tmp = M
tmp = tmp.subs(sub_list)
tmp /= newton*m
pprint(N(tmp,prec))
