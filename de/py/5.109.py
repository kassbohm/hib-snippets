r = 40*cm
(a, w) = (4 /s**2, 2 /s)
beta = 30*pi/180

(cb, sb) = (cos(beta), sin(beta))
dCB = Matrix([-r*cb, r*sb, 0])

vC = Matrix([r*w, 0, 0])
aC = Matrix([r*a, 0, 0])

vBx, vBy = var("vBx, vBy")
aBx, aBy = var("aBx, aBy")

W = Matrix([0,0,-w])
A = Matrix([0,0,-a])
vB = Matrix([vBx, vBy, 0])
aB = Matrix([aBx, aBy, 0])

# Find vBx and vBy:
eq = Eq( vB, vC + W.cross(dCB) )
sol = solve(eq, [vBx, vBy])
vBx, vBy = sol[vBx], sol[vBy]
vB = Matrix([vBx, vBy, 0])

pprint("\nvB / (m/s):")
tmp = vB
tmp /= m/s
tmp = iso_round(tmp,0.01)
pprint(tmp)

# Find aBx and aBy:
eq = Eq( aB, aC + A.cross(dCB) + W.cross(W.cross(dCB)) )
sol = solve(eq, [aBx, aBy])
aBx, aBy = sol[aBx], sol[aBy]
aB = Matrix([aBx, aBy, 0])

pprint("\naB / (m/s²):")
tmp = aB
tmp /= m/s**2
tmp = iso_round(tmp,0.01)
pprint(tmp)

# vB / (m/s):
# ⎡1.2 ⎤
# ⎢    ⎥
# ⎢0.69⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# aB / (m/s²):
# ⎡3.79⎤
# ⎢    ⎥
# ⎢0.59⎥
# ⎢    ⎥
# ⎣0.0 ⎦
