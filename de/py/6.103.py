from sympy.physics.units import *
from sympy import *

g, q, M, h, d = var("g, q, M, h, d")

sub_list = [
    ( g, S(981)/100 *m/s**2 ),
    ( q, 100 *newton/m  ),
    ( M, 600 *newton*m  ),
    ( h, 3 *m  ),
    ( d, 4 *m ),
]

half = S(1)/2
l = sqrt(h*h + d*d)

# masses:
m1g = q * h
m2g = q * l
m1 = m1g / g
m2 = m2g / g
theta1 = m1*h*h/12

# x1 = x_1''
x1, Ah, Av, C, Gh, Gv = var("x_1, A_h, A_v, C, G_h, G_v")

e1 = Eq( 2*theta1*x1/h,  M - half*h*Ah - half*h*Gh )
e2 = Eq( m1*x1,  Ah - Gh)
e3 = Eq( 0,  m1*g - Gv - Av)
e4 = Eq( 0,  h*Gh - d*Gv - d*C)
e5 = Eq( 2*m2*x1,  Gh )
e6 = Eq( 0, Gv + m2*g - C)

# solve linear system:
sol = solve([e1,e2,e3,e4,e5,e6], [x1, Ah, Av, C, Gh, Gv])
pprint("Gh, Gv, C / Newton:")
Gh, Gv, C = sol[Gh], sol[Gv], sol[C]
for s in [Gh, Gv, C]:
    tmp = s/newton
    tmp = tmp.subs(sub_list)
    pprint(N(tmp,3))
