# 9.22
b, c, d, l = var("b, c, d, l")

vAz, aAz = var("vAz, aAz")

sub_list = [
    ( b    ,   2  *m       ),
    ( c    ,   1  *m       ),
    ( d    ,   S(3)/2 *m   ),
    ( vAz  , - 3  *m/s     ),
    ( aAz  ,   0  *m/s**2  ),
]

l = sqrt(b*b + c*c)

dAB = Matrix([0, b, -c])

vA = Matrix([0, 0, vAz])
aA = Matrix([0, 0, aAz])

# ω: 3 unknowns:
wx, wy, wz = var("wx, wy, wz")
w = Matrix([wx, wy, wz])

# vBe: 1 unknown:
vBe = var("vBe")
ca = b/sqrt(b*b + d*d)
sa = d/sqrt(b*b + d*d)
e = Matrix([-sa, ca, 0])
vB = vBe * e

# equations to solve:
eq1 = Eq(w.dot(dAB), 0)
eq2 = Eq(vB, vA + w.cross(dAB))
sol = solve([eq1, eq2], [wx,wy,wz,vBe])
wx, wy, wz, vBe = sol[wx], sol[wy], sol[wz], sol[vBe]
w = Matrix([wx, wy, wz])
vB = vBe * e

pprint("\n\nω / (1/s):")
tmp = w
tmp =  tmp.subs(sub_list)
tmp = tmp / (1/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\n\nvBe / (m/s):")
tmp = vBe
tmp =  tmp.subs(sub_list)
tmp = tmp / (m/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\n\nvB / (m/s):")
tmp = vB
tmp =  tmp.subs(sub_list)
tmp = tmp / (m/s)
tmp = iso_round(tmp,0.001)
pprint(tmp)
         
# ω / (1/s):
# ⎡ 1.5 ⎤
# ⎢     ⎥
# ⎢0.225⎥
# ⎢     ⎥
# ⎣0.45 ⎦
#
#
# vBe / (m/s):
# 1.875
#
#
# vB / (m/s):
# ⎡-1.125⎤
# ⎢      ⎥
# ⎢ 1.5  ⎥
# ⎢      ⎥
# ⎣ 0.0  ⎦
