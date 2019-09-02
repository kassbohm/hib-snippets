from sympy.physics.units import *
from sympy import *

prec = 3

(lAB, rA, rZ) = (2*m, 1*m/2, 7*m/10)
(w, a) = (2 /s, 4 /s**2)
beta = asin( (rA + rZ)/lAB)
(cb, sb) = (cos(beta), sin(beta))

# Unknowns:
vAx, vAy, aAx, aAy = var("vAx, vAy, aAx, aAy")

# Angular velocity and acceleration:
# Body 1 = Wheel:
W1 = Matrix([0,0,w])
A1 = Matrix([0,0,a])

vA = Matrix([vAx, vAy, 0])
aA = Matrix([aAx, aAy, 0])
vO = Matrix([-rZ*w, 0, 0])
aO = Matrix([-rZ*a, 0, 0])

dOA = Matrix([0, rA, 0])

# Find vAx, vAy:
eq = Eq( vA, vO + W1.cross(dOA) )
sol = solve(eq, [vAx, vAy])
vAx, vAy = sol[vAx], sol[vAy]
vA = Matrix([vAx, vAy, 0])
pprint("\nvA / (m/s):")
tmp = vA
tmp /= m/s
tmp = N(tmp,prec)
pprint(tmp)

# Find aAx, aAy:
eq = Eq( aA, aO + A1.cross(dOA) + W1.cross(W1.cross(dOA)) )
sol = solve(eq, [aAx, aAy])
aAx, aAy = sol[aAx], sol[aAy]
aA = Matrix([aAx, aAy, 0])
pprint("\naA / (m/s²):")
tmp = aA
tmp /= m/s**2
tmp = N(tmp,prec)
pprint(tmp)

e = Matrix([-sb, cb, 0])

pprint("\nw2 / (1/s):")
vt = vA.dot(e)
w2 = vt/lAB
tmp = w2
tmp /= 1/s
pprint(N(tmp,prec))

# Body 2 = Rod BC:
# Angular velocity and acceleration:
W2 = Matrix([0,0,w2])
a2 = var("a2")
A2 = Matrix([0,0,a2])

# Change of coords:
R = Matrix([ [cb, sb, 0], [-sb, cb, 0], [0, 0, 1] ])
Rt = R.transpose()

ref = "(x\u0304, y\u0304, z\u0304)"
pprint("\n"+ref+"-comp's of the wheel-part.-velocity vA / (m/s):")
tmp = R*vA
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n"+ref+"-comp's of the vel. of a part. on the rod BC next to A / (m/s):")
vAb = Matrix([0, w2*lAB, 0])
tmp = vAb
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n"+ref+"comp's of the wheel-part.-velocity rel. to "+ref+" / (m/s):")
tmp = R*vA - vAb
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n(x,y,z)-comp's of the wheel-part.-velocity rel. to "+ref+" / (m/s):")
vxyz = Rt*(R*vA - vAb)
tmp = vxyz
tmp /= m/s
pprint(N(tmp,prec))

# unknown acceleration of the peg at A (on wheel) relative to (xbar, ybar, zbar)
ax = var("ax")
Axb = Matrix([ax, 0, 0])
Ax = Rt*Axb

dBA = Matrix([lAB*cb, lAB*sb, 0])
eq = Eq( aA, A2.cross(dBA) + W2.cross(W2.cross(dBA)) + 2*W2.cross(vxyz) + Ax)
sol = solve(eq, [ax, a2])
ax, a2 = sol[ax], sol[a2]

pprint("\na2 / (1/s²):")
tmp = a2
tmp /= 1/s**2
pprint(N(tmp,prec))

pprint("\nax / (m/s²):")
tmp = ax
tmp /= m/s**2
pprint(N(tmp,prec))
