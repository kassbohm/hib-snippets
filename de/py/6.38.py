from sympy.physics.units import *
from sympy import *

case = 1

prec = 3

newton = kg*m/s**2
kN = 1000*newton
g = S(981)/100 *m/s**2

mass, v, mu = var("mass, v, mu")
b, c, h = var("b, c, h")

sub_list = [
    (mass, 1500 *kg),
    (v, 80 * 1000*m / (3600 *s)),
    (mu, S(2)/10),
    (b, S(125)/100 *m),
    (c,  S(75)/100 *m),
    (h,  S(35)/100 *m),
    ]

a, NA, NB = var("a, NA, NB")
if case == 1:
    # For: mu NA = 0:
    eq1 = Eq(mass*a, mu*(0 + NB))
    eq2 = Eq(mass*g, NA + NB)
    eq3 = Eq(0, c*NB - b*NA - mu*h*(0 + NB))
elif case ==2:
    eq1 = Eq(mass*a, mu*(NA + NB))
    eq2 = Eq(mass*g, NA + NB)
    eq3 = Eq(0, c*NB - b*NA - mu*h*(NA + NB))

sol = solve([eq1, eq2, eq3], [a, NA, NB])

pprint("\nNA / kN:")
NA = sol[NA]
tmp = NA
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= kN
pprint(N(tmp,prec))

pprint("\nNB / kN:")
NB = sol[NB]
tmp = NB
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= kN
pprint(N(tmp,prec))

pprint("\na / (m/s²):")
a = sol[a]
tmp = a
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= m/s**2
pprint(N(tmp,prec))

pprint("\nt / s:")
t = v/a
tmp = t
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= s
pprint(N(tmp,prec))
