from sympy.physics.units import *
from sympy import *

g = 981*m/s**2/100
newton = kg*m/s**2

(r, R) = ( 0.1 *m, 0.2 *m)
G = 300 *newton
(mass, Mass) = (10 *kg, G / g)
J = (0.13 *m)**2 * Mass

# (X, x, p) = (X'', x'', phi'')
X, x, p, H, S = var("X, x, p, H, S ")

e1 = Eq(X, r*p)
e2 = Eq(x, (r+R)*p)
e3 = Eq(mass*x, mass*g - S)
e4 = Eq(Mass*X, S - H)
e5 = Eq(J*p, R*S + r*H)

sol = solve([e1,e2,e3,e4,e5], [X, x, p, H, S])
x = sol[x]

pprint("\nx''/ (m/s²):")
tmp = x
tmp /= m/s**2
pprint(N(tmp,3))

pprint("\nx'(2 s) / (m/s):")
v = x* 2*s
tmp = v
tmp /= m/s
pprint(N(tmp,3))

# (X, x, p) = (X', x', phi')
X, x, p, H, S = var("X, x, p, H, S ")

t = 2*s
e1 = Eq(X, r*p)
e2 = Eq(x, (r+R)*p)
e3 = Eq(mass*x, (mass*g - S)*t)
e4 = Eq(Mass*X, (S - H)*t)
e5 = Eq(J*p, (R*S + r*H)*t)

sol = solve([e1,e2,e3,e4,e5], [X, x, p, H, S])

x = sol[x]
pprint("\nx'(2 s) / (m/s):")
tmp = x
tmp /= m/s
pprint(N(tmp,3))
