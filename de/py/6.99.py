from sympy.physics.units import *
from sympy import *

prec = 3

g = var("g")
M, l = var("M, l")
beta = var("beta")

sub_list = [
    ( g, 9.81  *m/s/s ),
    ( M, 2  *kg ),
    ( l, S(3)/10  *m ),
    ( beta, 30*pi/180),
]

I = S(1)/12 *M*l*l
cb, sb = cos(beta), sin(beta)

ax,ay, aBx,aBy,aB, S,alpha = var("ax,ay, aBx,aBy,aB, S,alpha")
unks = [ax,ay, aBx,aBy,aB, S,alpha]

e1 = Eq( M*ax , S*cb )
e2 = Eq( M*ay , -M*g + S*sb )
e3 = Eq( I*alpha , l/2*S*sb )
e4 = Eq( aBx, -aB*sb )
e5 = Eq( aBy,  aB*cb )
e6 = Eq( aBx,  ax )
e7 = Eq( aBy,  ay + alpha*l/2 )

eqs = [e1,e2,e3,e4,e5,e6,e7]

# solve linear system:
sol = solve(eqs, unks)
ax, ay = sol[ax], sol[ay]

pprint("\nax / (m/s²):")
tmp = ax
tmp = tmp.subs(sub_list)
tmp /= m/s**2
pprint(N(tmp,prec))

pprint("\nay / (m/s²):")
tmp = ay
tmp = tmp.subs(sub_list)
tmp /= m/s**2
pprint(N(tmp,prec))
