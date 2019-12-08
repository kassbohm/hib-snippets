prec = 3

rB, l = 150 *mm, 500 *mm
w, a  = 8 /s, 16 /s**2
phi   = 30 *deg
w = Matrix([0,0,-w])
a = Matrix([0,0,-a])


cp, sp = cos(phi), sin(phi)
er = Matrix([ -cp, sp, 0 ])
el = Matrix([  sp, cp, 0 ])
dBA = l*el

vAx = var("vAx")
vA = Matrix([vAx, 0, 0])

vBx, vBy = var("vBx, vBy")
vB = Matrix([vBx, vBy, 0])

wS = var("wS")
WS = Matrix([0,0,-wS])

# ---

eq = Eq(vB, w.cross(rB*er))
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
WS = Matrix([0,0,-wS])

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
