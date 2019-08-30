from sympy.physics.units import *
from sympy import *

r = S(3)/10 *meter
v = 2 * meter/second
a = 3 * meter/second/second

c1, s1 = cos(50 *degrees), sin(50 *degrees)

omega = v/(r*s1)
alpha = 1/s1*(a/r - omega**2*c1)

pprint("\nomega:")
tmp = omega
pprint(N(tmp,3))

pprint("\nalpha:")
tmp = alpha
pprint(N(tmp,3))

# omega:
# 8.7
# ───
#  s
#
# alpha:
# -50.5
# ──────
#    2
#   s
