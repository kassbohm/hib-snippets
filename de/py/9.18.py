from sympy.physics.units import *
from sympy import *

# Given quantities:
b, c, d, h, wAC, aAC = var("b, c, d, h, ωAC, αAC")

sub_list = [
    (b, 2 *cm),
    (c, 3 *cm),
    (d, S(3)/2 *cm),
    (h, 6 *cm),
    (wAC, 8 /s),
    (aAC, 0),         # 9.18
    # (aAC, 6 /s**2),     # 9.19
]

# AC:
w1 = Matrix([0, 0, wAC])
a1 = Matrix([0, 0, aAC])
# BD:
w2x, a2x = var("ω2x, α2x")
w2 = Matrix([w2x, 0, 0])
a2 = Matrix([a2x, 0, 0])
# AB:
w3x, w3y, w3z = var("ω3x, ω3y, ω3z")
a3x, a3y, a3z = var("α3x, α3y, α3z")
w3 = Matrix([w3x, w3y, w3z])
a3 = Matrix([a3x, a3y, a3z])

# A:
vAx,vAy,vAz, aAx,aAy,aAz = var("vAx,vAy,vAz, aAx,aAy,aAz")
vA = Matrix([vAx, vAy, vAz])
aA = Matrix([aAx, aAy, aAz])
# B:
vBx,vBy,vBz, aBx,aBy,aBz = var("vBx,vBy,vBz, aBx,aBy,aBz")
vB = Matrix([vBx, vBy, vBz])
aB = Matrix([aBx, aBy, aBz])

# --- a ---
dCA = Matrix([d, 0, 0])
dDB = Matrix([0, b, 0])
dBA = Matrix([-c, -b, h])

# --- b ---
# 6 equations for 6 unknowns:
unk = [vAx,vAy,vAz, aAx,aAy,aAz]
eq1 = Eq(vA, w1.cross(dCA))
eq2 = Eq(aA, a1.cross(dCA) + w1.cross(w1.cross(dCA)))
sol = solve([eq1, eq2], unk)
vA = Matrix([sol[vAx], sol[vAy], sol[vAz]])
aA = Matrix([sol[aAx], sol[aAy], sol[aAz]])

# --- c ---
# 6 equations for 8 unknowns:
unk = [vBx,vBy,vBz,aBx,aBy,aBz]
eq1 = Eq(vB, w2.cross(dDB))
eq2 = Eq(aB, a2.cross(dDB) + w2.cross(w2.cross(dDB)))
sol = solve([eq1, eq2], unk, dict=True)
sol = sol[0]
vB = Matrix([sol[vBx], sol[vBy], sol[vBz]])
aB = Matrix([sol[aBx], sol[aBy], sol[aBz]])

# --- d ---
# Velocities and Angular Velocities:
unk = [w2x, w3x,w3y,w3z]
eq1 = Eq(vA, vB + w3.cross(dBA))
eq2 = Eq(w3.dot(dBA),0)
sol = solve([eq1,eq2], unk, dict=True)
sol = sol[0]
w3 = Matrix([sol[w3x], sol[w3y], sol[w3z]])
w2x_sol = sol[w2x]
w2 = Matrix([w2x_sol, 0, 0])

# --- e ---
# Accelerations and Angular Accelerations:
eq1 = Eq(aA, aB + a3.cross(dBA) + w3.cross(w3.cross(dBA)))
eq2 = Eq(a3.dot(dBA),0)
eqns = [eq1, eq2]
unks = [a2x,a3x,a3y,a3z]
sol = solve(eqns, unks, dict=True)
sol = sol[0]
a2x,a3x,a3y,a3z = sol[a2x],sol[a3x],sol[a3y],sol[a3z]
a2x = a2x.subs(w2x, w2x_sol).simplify()
a3x = a3x.subs(w2x, w2x_sol).simplify()
a3y = a3y.subs(w2x, w2x_sol).simplify()
a3z = a3z.subs(w2x, w2x_sol).simplify()
a3 = Matrix([a3x, a3y, a3z])

pprint("\n\nω3:")
tmp = w3
pprint(tmp)
pprint("\n\nω3 / (1/s):")
tmp = w3.subs(sub_list)
tmp /= 1/s
pprint(tmp)

pprint("\n\na3:")
tmp = a3
pprint(tmp)
pprint("\n\na3 / (1/s²):")
tmp = a3.subs(sub_list)
tmp /= 1/s**2
pprint(tmp)
