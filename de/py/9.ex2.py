from sympy.physics.units import *
from sympy import *

(wS, wSp) = ( 10 /s, 6 /s**2 )
(wN, wNp) = (  3 /s, 2 /s**2 )
(wP, wPp) = (  5 /s, 4 /s**2 )

theta = 60 *pi/180
(ct, st) = (cos(theta), sin(theta))

w = Matrix([-wN, wS*st, wP + wS*ct])

pprint("\nω / (1/s):")
tmp = w / (1/s)
pprint(N(tmp,3))

wp  = Matrix([-wNp, wSp*st, wSp*ct + wPp])
O = Matrix([-wN, 0, wP])
wp += O.cross(w) + O.cross(Matrix([-wN, 0, 0]))

pprint("\nω' / (1/s²):")
tmp = wp / (1/s/s)
pprint(N(tmp,3))

# ω / (1/s):
# ⎡-3.0⎤
# ⎢    ⎥
# ⎢8.66⎥
# ⎢    ⎥
# ⎣10.0⎦
#
# ω' / (1/s²):
# ⎡-45.3⎤
# ⎢     ⎥
# ⎢ 5.2 ⎥
# ⎢     ⎥
# ⎣-19.0⎦
