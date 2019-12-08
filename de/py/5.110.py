rB, l = 150 *mm, 500 *mm
w, a  = 8 /s, 16 /s**2
phi   = 30 *deg
w = Matrix([0,0,-w])
a = Matrix([0,0,-a])

cp, sp = cos(phi), sin(phi)
er = Matrix([ -cp, sp, 0 ])
el = Matrix([  sp, cp, 0 ])
dOB = rB*er
dBA =  l*el

# ---

vAx, aAx = var("vAx, aAx")
vA = Matrix([vAx, 0, 0])
aA = Matrix([aAx, 0, 0])

vBx, vBy = var("vBx, vBy")
aBx, aBy = var("aBx, aBy")
vB = Matrix([vBx, vBy, 0])
aB = Matrix([aBx, aBy, 0])

wS, aS = var("wS, aS")
WS = Matrix([0,0,wS])
AS = Matrix([0,0,aS])

# ---

eq = Eq(vB, w.cross(dOB))
sol = solve(eq, [vBx, vBy])
vBx, vBy = sol[vBx], sol[vBy]
vB = Matrix([vBx, vBy, 0])

pprint("\nvB / (m/s):")
tmp = vB
tmp /= (m/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)

eq = Eq(vA, vB + WS.cross(dBA))
sol = solve(eq, [vAx, wS])
vAx, wS = sol[vAx], sol[wS]
vA = Matrix([vAx, 0, 0])
WS = Matrix([0,0,wS])

pprint("\nvA / (m/s):")
tmp = vA
tmp /= (m/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)

pprint("\nWS / (1/s):")
tmp = WS
tmp /= (1/s)
tmp = iso_round(tmp,0.01)
pprint(tmp)

eq = Eq( aB, a.cross(dOB) + w.cross(w.cross(dOB)) )
sol = solve(eq, [aBx, aBy])
aBx, aBy = sol[aBx], sol[aBy]
aB = Matrix([aBx, aBy, 0])

pprint("\naB / (m/s²):")
tmp = aB
tmp /= (m/s**2)
tmp = iso_round(tmp,0.01)
pprint(tmp)

eq = Eq( aA, aB + AS.cross(dBA) + WS.cross(WS.cross(dBA)) )
sol = solve(eq, [aAx, aS])
aAx, aS = sol[aAx], sol[aS]
aA = Matrix([aAx, 0, 0])
AS = Matrix([0,0,aS])

pprint("\naA / (m/s²):")
tmp = aA
tmp /= (m/s**2)
tmp = iso_round(tmp,0.01)
pprint(tmp)

pprint("\nAS / (1/s²):")
tmp = AS
tmp /= (1/s**2)
tmp = iso_round(tmp,0.01)
pprint(tmp)

# vB / (m/s):
# ⎡0.6 ⎤
# ⎢    ⎥
# ⎢1.04⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# vA / (m/s):
# ⎡2.4⎤
# ⎢   ⎥
# ⎢0.0⎥
# ⎢   ⎥
# ⎣0.0⎦
#
# WS / (1/s):
# ⎡ 0.0 ⎤
# ⎢     ⎥
# ⎢ 0.0 ⎥
# ⎢     ⎥
# ⎣-4.16⎦
#
# aB / (m/s²):
# ⎡9.51 ⎤
# ⎢     ⎥
# ⎢-2.72⎥
# ⎢     ⎥
# ⎣ 0.0 ⎦
#
# aA / (m/s²):
# ⎡-12.48⎤
# ⎢      ⎥
# ⎢ 0.0  ⎥
# ⎢      ⎥
# ⎣ 0.0  ⎦
#
# AS / (1/s²):
# ⎡ 0.0 ⎤
# ⎢     ⎥
# ⎢ 0.0 ⎥
# ⎢     ⎥
# ⎣40.82⎦
