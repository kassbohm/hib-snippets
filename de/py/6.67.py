from sympy.physics.units import *
from sympy import *

prec = 3
newton = kg*m/s**2
g = S(981)/100 *m/s**2

mass, l = var("mass, l")

sub_list = [
    (mass, 4 *kg),
    (l, 2 *m),
    ]

theta = mass*l**2/12

a, alpha = var("a, alpha")
eq1 = Eq( mass*a, mass*g - mass*g / 2 )
eq2 = Eq( theta*alpha, l*mass*g/4 )

sol = solve([eq1, eq2], [a, alpha])
a, alpha = sol[a], sol[alpha]

pprint("\na / (m/s²):")
tmp = a
tmp = tmp.subs(sub_list)
tmp /= m/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nalpha / (1/s²):")
tmp = alpha
tmp = tmp.subs(sub_list)
tmp /= 1/s**2
tmp = N(tmp,prec)
pprint(tmp)
