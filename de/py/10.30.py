from sympy.physics.units import *
from sympy import *

M = var("mass")
g = 9.81 *m/s/s
(d , h) = (1.5 *m , 2.5 *m)

r = d/2

pprint("\nv1 / (m/s):")
v1 = sqrt(2*g*h)
tmp = v1
tmp /= m/s
pprint(N(tmp,3))

J_S = M*r**2 / 4
J_A = J_S + M*r*r

v2 = var("v2")

eq = Eq(r*M*v1, J_A*v2/r)

sol = solve(eq,v2)[0]
pprint("\nv2 / (m/s):")
tmp = sol
tmp /= m/s
pprint(N(tmp,3))

# v1 / (m/s):
# 7.00
#
# v2 / (m/s):
# 5.60
