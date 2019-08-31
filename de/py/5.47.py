from sympy.physics.units import *
from sympy import *

t = var("t")
ω, α = var("ω, α")
C, S = var("C, S")
r, l = 3 *m, 5 *m

theta = Function("theta")(t)
s = sqrt(r*r + l*l - 2*r*l*cos(theta))
s = s.simplify()

v = s.diff(t)
a = v.diff(t)

pprint("\nv:")
tmp = v.subs(theta.diff(t),ω)
tmp = tmp.subs(sin(theta),S)
tmp = tmp.subs(cos(theta),C)
pprint(tmp)


pprint("\na:")
tmp = a.subs(theta.diff(t),ω)
tmp = tmp.subs(sin(theta),S)
tmp = tmp.subs(cos(theta),C)
pprint(tmp)
