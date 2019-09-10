from sympy.physics.units import *
from sympy import *

prec = 3
cm = m/100

l = 10*cm
(vA, aA) = (5 *cm/s, 7 *cm/s**2)
beta = 60*pi/180
vA = Matrix([0, -vA, 0])
aA = Matrix([0, -aA, 0])

# ---

(cb, sb) = (cos(beta), sin(beta))
dAB = Matrix([l*cb, -l*sb, 0])

w, a = var("w, a")
W = Matrix([0,0,w])
A = Matrix([0,0,a])
vBx, aBx = var("vBx, aBx")
vB = Matrix([vBx, 0, 0])
aB = Matrix([aBx, 0, 0])

# Find vBx and w:
eq = Eq( vB, vA + W.cross(dAB) )
sol = solve(eq, [vBx, w])
vBx, w = sol[vBx], sol[w]
vB = Matrix([vBx, 0, 0])
W = Matrix([0,0,w])

pprint("\nvB / (cm/s):")
tmp = vB
tmp /= cm/s
tmp = N(tmp,prec)
pprint(tmp)

pprint("\nw / (1/s):")
tmp = W
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)

# Find aBx and alpha:
eq = Eq( aB, aA + A.cross(dAB) + W.cross(W.cross(dAB)) )
sol = solve(eq, [aBx, a])
aBx, a = sol[aBx], sol[a]
aB = Matrix([aBx, 0, 0])
A = Matrix([0,0,a])

pprint("\naB / (cm/s²):")
tmp = aB
tmp /= cm/s**2
tmp = N(tmp,prec)
pprint(tmp)

pprint("\na / (1/s²):")
tmp = A
tmp /= (1/s**2)
tmp = N(tmp,prec)
pprint(tmp)
            
# vB / (cm/s):
# ⎡8.66⎤
# ⎢    ⎥
# ⎢ 0  ⎥
# ⎢    ⎥
# ⎣ 0  ⎦
#
# w / (1/s):
# ⎡  0   ⎤
# ⎢      ⎥
# ⎢  0   ⎥
# ⎢      ⎥
# ⎣second⎦
#
# aB / (cm/s²):
# ⎡-7.88⎤
# ⎢     ⎥
# ⎢  0  ⎥
# ⎢     ⎥
# ⎣  0  ⎦
#
# a / (1/s²):
# ⎡  0   ⎤
# ⎢      ⎥
# ⎢  0   ⎥
# ⎢      ⎥
# ⎣-0.332⎦
