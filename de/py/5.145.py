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
pprint("\n\n--- Solution A ---")
# Angular velocity and acceleration:
W2 = Matrix([0,0,w2])
a2 = var("a2")
A2 = Matrix([0,0,a2])

# Change of coords:
R = Matrix([ [cb, sb, 0], [-sb, cb, 0], [0, 0, 1] ])
Rt = R.transpose()

pprint("\nAll velocities in (m/s):")

ref = "(x\u0304,y\u0304,z\u0304)"
pprint("\n"+ref+"-comp's of the velocity of A relative to the ground:")
tmp = R*vA
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n"+ref+"-comp's of the vel. of a part. on BC next to A rel. to ground:")
vAb = Matrix([0, w2*lAB, 0])
tmp = vAb
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n"+ref+"comp's of vAB, the velocity of A relative to rotating "+ref+":")
tmp = R*vA - vAb
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n(x,y,z)-comp's of vAB:")
vAB = Rt*(R*vA - vAb)
tmp = vAB
tmp /= m/s
pprint(N(tmp,prec))

# unknown acceleration of the peg at A (on wheel) relative to (xbar, ybar, zbar)
a = var("a")
Axb = Matrix([a, 0, 0])
aAB = Rt*Axb

dBA = Matrix([lAB*cb, lAB*sb, 0])
eq = Eq( aA, A2.cross(dBA) + W2.cross(W2.cross(dBA)) + 2*W2.cross(vAB) + aAB)
sol = solve(eq, [a, a2])
a, a2 = sol[a], sol[a2]

pprint("\nα2 / (1/s²):")
tmp = a2
tmp /= 1/s**2
pprint(N(tmp,prec))

pprint("\na / (m/s²):")
tmp = a
tmp /= m/s**2
pprint(N(tmp,prec))

pprint("\n\n--- Solution B ---")
a , a2= var("a, a2")
A2 = Matrix([0,0,a2])

pprint("\n(x,y,z)-comp's of v_AB / (m/s):")
vAB = vA - W2.cross(dBA)
tmp = vAB
tmp /= m/s
pprint(N(tmp,prec))

pprint("\n(x,y,z)-comp's of a_AB:")
aAB = Matrix([a*cb, a*sb, 0])
tmp = aAB
pprint(N(tmp,prec))

eq = Eq( aA, A2.cross(dBA) + W2.cross(W2.cross(dBA)) + 2*W2.cross(vAB) + aAB)
sol = solve(eq, [a, a2])
a, a2 = sol[a], sol[a2]

pprint("\nα2 / (1/s²):")
tmp = a2
tmp /= 1/s**2
pprint(N(tmp,prec))

pprint("\na / (m/s²):")
tmp = a
tmp /= m/s**2
pprint(N(tmp,prec))

# vA / (m/s):
# ⎡-2.4⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# aA / (m/s²):
# ⎡-4.8⎤
# ⎢    ⎥
# ⎢-2.0⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# w2 / (1/s):
# 0.720
#
#
# --- Solution A ---
#
# All velocities in (m/s):
#
# (x̄,ȳ,z̄)-comp's of the velocity of A relative to the ground:
# ⎡-1.92⎤
# ⎢     ⎥
# ⎢1.44 ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# (x̄,ȳ,z̄)-comp's of the vel. of a part. on BC next to A rel. to ground:
# ⎡ 0  ⎤
# ⎢    ⎥
# ⎢1.44⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# (x̄,ȳ,z̄)comp's of vAB, the velocity of A relative to rotating (x̄,ȳ,z̄):
# ⎡-1.92⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# (x,y,z)-comp's of vAB:
# ⎡-1.54⎤
# ⎢     ⎥
# ⎢-1.15⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# α2 / (1/s²):
# 2.02
#
# a / (m/s²):
# -4.00
#
#
# --- Solution B ---
#
# (x,y,z)-comp's of v_AB / (m/s):
# ⎡-1.54⎤
# ⎢     ⎥
# ⎢-1.15⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# (x,y,z)-comp's of a_AB:
# ⎡0.8⋅a⎤
# ⎢     ⎥
# ⎢0.6⋅a⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# α2 / (1/s²):
# 2.02
#
# a / (m/s²):
# -4.00
