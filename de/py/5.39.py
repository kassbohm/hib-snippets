from sympy.physics.units import *
from sympy import *

prec = 4

r = 3*m/10
v = 2 * m/s
a = 3 * m/s**2

p1 = rad(50)
p1 = 50*pi/180
c1 = cos(p1)
s1 = sin(p1)

omega = v/(r*s1)
pprint("\nomega / (1/s):")
tmp = omega
tmp /= 1/s
pprint(N(tmp,prec))


pprint("\nalpha / (1/s²):")
alpha  = N(a/r)
alpha -= N(omega**2*c1)
alpha /= s1
tmp = alpha
tmp /= (1/s**2)
pprint(N(tmp,prec))
