la, lb, ha = var("la, lb, ha")
alpha = var("alpha")
rD, v = var("rD, v")

sub_list=[
    ( la,     2*m           ),
    ( lb,     1*m           ),
    ( ha,     3*m/2         ),
    ( alpha, 45 *pi/180     ),
    ( rD,     5*m  /100     ),
    ( v,     25*m/s/100     ),
    ]

wD = v/rD

# to reproduce Hibbeler's (wrong) solution:
# wD *= 2

ca, sa = cos(alpha), sin(alpha)
t1 = Matrix([0, 0, ha])
t2 = la * Matrix([ca, 0,  sa])
t3 = lb * Matrix([sa, 0, -ca])
rC = t1 + t2 + t3
WD = Matrix([0,0,wD])

vC = WD.cross(rC)
aC = WD.cross(WD.cross(rC))

pprint("\nvC / (m/s):   ")
tmp = vC
tmp = tmp.subs(sub_list)
tmp /= m/s
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\naC / (m/s):   ")
tmp = aC
tmp = tmp.subs(sub_list)
tmp /= m/s**2
tmp = iso_round(tmp,0.1)
pprint(tmp)

# vC / (m/s):
# ⎡0.0 ⎤
# ⎢    ⎥
# ⎢10.6⎥
# ⎢    ⎥
# ⎣0.0 ⎦
#
# aC / (m/s):
# ⎡-53.0⎤
# ⎢     ⎥
# ⎢ 0.0 ⎥
# ⎢     ⎥
# ⎣ 0.0 ⎦
