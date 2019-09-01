from sympy.physics.units import *
from sympy import *

newton = kg*m/s**2

M, w = var("M, omega")
h, l, theta = var("h, l, theta")
g = var("g")

h2 = h*h
(Jx, Jy, Jz) = (M*h*h/12, M*h*h/12, M*h*h/6)

# unknowns:
Ax, Ay, Bx, By = var("Ax, Ay, Bx, By")

sub_list=[
    (theta, 30*pi/180),
    (l, 3 *m),
    (h, S(1)/2 *m),
    (M, 20 *kg),
    (w, 25 /s),
    (g, S(981)/100 *m/s**2)
    ]

# Shortcuts:
s2theta = sin(2*theta)
w2 = w*w

eq1 = Eq(w2 * s2theta * (Jy - Jz)/2, -l/2*By + l/2*Ay)
eq2 = Eq(Ax, Bx)
eq3 = Eq(0, Bx+Ax)
eq4 = Eq(M*g, By+Ay)

eqns = [eq1, eq2, eq3, eq4]

pprint("\n --- Equations: ---")
for e in eqns:
    pprint("\n")
    pprint(e)

unknowns = [Ax, Ay, Bx, By]
sol = solve(eqns, unknowns)

Ax,Ay, Bx,By = sol[Ax],sol[Ay], sol[Bx],sol[By]

pprint("\n --- Reactions: ---")
pprint("\n --- rhs in the book, i.e. lhs here: ---")
pprint("\nBy:\n\n")
pprint(By)
pprint("\nBy / N: \n")
tmp = By.simplify()
tmp = By.subs(sub_list)
tmp /= newton
pprint(N(tmp,3))

pprint("\n --- lhs in the book, i.e. rhs here: ---")
pprint("\nAy:\n\n")
pprint(Ay)
pprint("\nAy / N:")
tmp = Ay
tmp = tmp.subs(sub_list)

pprint(N(tmp,3))
