from sympy.physics.units import *
from sympy import *

# Given:
a = var("a")
w1, a1 = var("w1, a1")

sub_list = [
    (a, 4*m/10),
    (w1, 3/s),
    (a1, 4/s**2),
    ]

# Unknowns:
w2, a2 = var("w2, a2")
vC, aC = var("vC, aC")
prec = 3

# Omega:
W = Matrix([0, 0, -w1])
O = Matrix([0, 0, -w2])
Wp = Matrix([0, 0, -a1])
Op = Matrix([0, 0, -a2])

# Velocity and acceleration
# of a particle on rigid body 1
# at C and relative to (x,y,z):
VC = Matrix([vC, 0, 0])
AC = Matrix([aC, 0, 0])

dAC = Matrix([a, a, 0])
dDC = Matrix([a, 0, 0])

# Equations to solve:
lhs = W.cross(dAC)
rhs = O.cross(dDC) + VC
eq = Eq(lhs,rhs)
sol = solve(eq, [w2, vC])

pprint("\nw2 / (1/s):")
w2_sol = sol[w2]
tmp = w2_sol
tmp = tmp.subs(sub_list)
tmp /= 1/s
pprint(tmp)

pprint("\nvC / (m/s):")
vC_sol = sol[vC]
tmp = vC_sol
tmp = tmp.subs(sub_list)
tmp /= m/s
pprint(tmp)

lhs = Wp.cross(dAC) - w1**2*dAC
rhs = Op.cross(dDC) + O.cross(O.cross(dDC)) + 2*O.cross(VC) + AC
eq = Eq(lhs,rhs)
sol = solve(eq, [a2, aC])

pprint("\na2 / (1/s²):")
tmp = sol[a2]
tmp = tmp.subs(w2, w2_sol)
tmp = tmp.subs(vC, vC_sol)
tmp = tmp.subs(sub_list)
tmp /= 1/s**2
pprint(tmp)

pprint("\naC / (m/s²):")
tmp = sol[aC]
tmp = tmp.subs(w2, w2_sol)
tmp = tmp.subs(vC, vC_sol)
tmp = tmp.subs(sub_list)
tmp /= m/s**2
pprint(tmp)

# w2 / (1/s):
# 3
#
# vC / (m/s):
# 6/5
#
# a2 / (1/s²):
# -5
#
# aC / (m/s²):
# 8/5
