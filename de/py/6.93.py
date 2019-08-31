from sympy.physics.units import *
from sympy import *

prec = 3

kilo = 1000
newton = kg*m/s**2
kN = kilo*newton

M, k, a, ra, ri, muH, muG = var("M, k, a, ra, ri, muH, muG")

sub_list = [
    ( M , 500  *kg ),
    ( k , S(13)/10  *m ),
    ( a , 1  *m/s/s ),
    ( ra , S(16)/10  *m ),
    ( ri , S(8)/10  *m ),
    ( muH , S(5)/10 ),
    ( muG , S(4)/10 ),
    ]

# grav. acceleration:
g = 10 *m/s/s
theta = k*k*M

# Unknowns:
H, S, x, phi = var("H, S, x, phi")

# equilibrium conditions:
e1 = Eq( theta*phi , ra*H - ri*S )
e2 = Eq( M*x , H - S )
e3 = Eq( x + ra*phi , a )
e4 = Eq( x , -ri*phi )

# solve linear system:
sol = solve([e1,e2,e3,e4], [H, S, x, phi])

S = sol[S]
H = sol[H]
phi = sol[phi]

pprint("\nt=0:")

pprint("\nSeilkraft S / kN:")
tmp = S
tmp = tmp.subs(sub_list)
tmp /= kN
pprint(N(tmp,prec))

pprint("\nphi'' / (1/s²):")
tmp = phi
tmp = tmp.subs(sub_list)
tmp /= 1/s**2
pprint(N(tmp,prec))

pprint("\nHaftkraft H / kN:")
tmp = H
tmp = tmp.subs(sub_list)
tmp /= kN
pprint(N(tmp,prec))

pprint("\nNormalkraft N / kN:")
tmp = M*g
tmp = tmp.subs(sub_list)
tmp /= kN
pprint(N(tmp,prec))
