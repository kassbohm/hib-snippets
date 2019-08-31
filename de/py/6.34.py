from sympy.physics.units import *
from sympy import *

prec = 5

newton = kg*m/s**2
g = S(981)/100 *m/s**2

mass, a, mu = var("mass, a, mu")

sub_list = [
    (mass, 800 *kg),
    (a, S(5)/10 *m/s**2),
    (mu, S(1)/10),
    ]

G = mass*g
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
tmp /= newton
pprint(N(tmp,prec))

pprint("\nS / Newton:")
S = sol[S]
tmp = S
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= newton
pprint(N(tmp,prec))

pprint("\nphi / deg:")
phi = asin(mu * Nv / S)
tmp = phi
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp = deg(tmp)
pprint(N(tmp,prec))
pprint("\ntheta / deg:")
tmp = 45 - tmp
pprint(N(tmp,prec))
