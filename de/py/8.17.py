from sympy.physics.units import *
from sympy import *

prec = 3

mm = m/1000

g = 981*m/s**2/100
mass = 70 *kg
r = 300 *mm
k0 = 125 *mm
alpha = 30 *pi/180

(muh, mug) = ( S(4)/10, S(3)/10 )
t = 2 *second

sa, ca, ta = sin(alpha), cos(alpha), tan(alpha)

J0 = mass *k0*k0
JA = J0 + mass*r*r

# Angular acceleration:
p = mass *g*sa*r / JA
# Angular velocity:
w = p * t

mu = (1 - mass*r*r/JA) *ta

pprint("\nμ:")
tmp = mu
pprint(N(tmp,prec))

pprint("\nω / (1/s):")
tmp = w
tmp /= 1/s
pprint(N(tmp,prec))
