from sympy.physics.units import *
from sympy import *
from sympy import N as Num

prec = 3
newton = kg*m/s**2
mm = m / 1000
g = S(981)/100 *m/s**2

mass, k, F, mu, r, h = var("mass, k, F, mu, r, h")

phi = atan(h/r)
theta = k**2*mass

sub_list = [
    (mass, 20 *kg),
    (k, 90 *mm),
    (F, 30 *newton),
    (mu, S(2)/10),
    (r, 125 *mm),
    (h, 300 *mm),
    ]

pprint("\nphi / °:")
tmp = phi
tmp = tmp.subs(sub_list)
tmp = tmp *180/pi
pprint(N(tmp,prec))

cc, ss = cos(phi), sin(phi)

N, S, a = var("N, S, a")
eq1 = Eq( N, S*cc )
eq2 = Eq( mu*N + mass*g + F, S*ss )
eq3 = Eq( theta*a, r*F - r*mu*N )

sol = solve([eq1, eq2, eq3], [N, S, a])
N, S, a = sol[N], sol[S], sol[a]

pprint("\nN / Newton:")
tmp = N
tmp = tmp.subs(sub_list)
tmp /= newton
tmp = Num(tmp,prec)
pprint(tmp)

pprint("\nS / Newton:")
tmp = S
tmp = tmp.subs(sub_list)
tmp /= newton
tmp = Num(tmp,prec)
pprint(tmp)

pprint("\na / (1/s²):")
tmp = a
tmp = tmp.subs(sub_list)
tmp /= 1/s**2
tmp = Num(tmp,prec)
pprint(tmp)
