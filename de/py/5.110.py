prec = 3

(r, rA, l) = (3*cm, 2*cm, 8*cm)
(w, a) = (6 /s, 12 /s**2)
beta = 60*pi/180
(cb, sb) = (cos(beta), sin(beta))

w1 = Matrix([0,0,w])
a1 = Matrix([0,0,a])

vAx, vAy, aAx, aAy = var("vAx, vAy, aAx, aAy")
vBx, aBx = var("vBx, aBx")

vA = Matrix([vAx, vAy, 0])
aA = Matrix([aAx, aAy, 0])

vB = Matrix([vBx, 0, 0])
aB = Matrix([aBx, 0, 0])

aP = Matrix([0, r*w**2, 0])

h = r - rA
dPA = Matrix([0, h, 0])
dAB = Matrix([l*cb, l*sb, 0])

# Rigid body 1:
# Find vAx and vAy:
eq = Eq( vA, w1.cross(dPA) )
sol = solve(eq, [vAx, vAy])
vAx, vAy = sol[vAx], sol[vAy]
vA = Matrix([vAx, vAy, 0])

pprint("\nvA / (cm/s):")
tmp = vA
tmp /= cm/s
tmp = N(tmp,prec)
pprint(tmp)

# Find aAx and aAy:
eq = Eq( aA, aP + a1.cross(dPA) + w1.cross(vA) )
sol = solve(eq, [aAx, aAy])
aAx, aAy = sol[aAx], sol[aAy]
aA = Matrix([aAx, aAy, 0])

pprint("\naA / (cm/s²):")
tmp = aA
tmp /= cm/s**2
tmp = N(tmp,prec)
pprint(tmp)

# Rigid body 2:
w2, a2 = var("w2, a2")
W2 = Matrix([0,0,w2])
A2 = Matrix([0,0,a2])

# Find vBx and w2:
eq = Eq( vB, vA + W2.cross(dAB) )
sol = solve(eq, [vBx, w2])
vBx, w2 = sol[vBx], sol[w2]
vB = Matrix([vBx, 0, 0])

pprint("\nvB / (cm/s):")
tmp = vB
tmp /= cm/s
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nw2 / (1/s):")
W2 = Matrix([0,0,w2])
tmp = W2
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)

# Find aBx and a2:
eq = Eq( aB, aA + A2.cross(dAB) + W2.cross(W2.cross(dAB)) )
sol = solve(eq, [aBx, a2])
aBx, a2 = sol[aBx], sol[a2]

pprint("\naB / (cm/s²):")
aB = Matrix([aBx, 0, 0])
tmp = aB
tmp /= cm/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\na2 / (1/s²):")
A2 = Matrix([0,0,a2])
tmp = A2
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)

# vA / (cm/s):
# ⎡-6.0⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# aA / (cm/s²):
# ⎡-12.0⎤
# ⎢     ⎥
# ⎢72.0 ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# vB / (cm/s):
# ⎡-6.0⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# w2 / (1/s):
# ⎡0⎤
# ⎢ ⎥
# ⎢0⎥
# ⎢ ⎥
# ⎣0⎦
#
# aB / (cm/s²):
# ⎡113.0⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# a2 / (1/s²):
# ⎡  0  ⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣-18.0⎦
