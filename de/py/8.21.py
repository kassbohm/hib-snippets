from sympy.physics.units import *
from sympy import *

prec=3

newton = kg*m/s**2
mm = m/1000

mass = 12 *kg
r = 200 *mm
w0 = 20 /s
mu = S(4)/10

J = mass*r*r / 2
k = 2 *mu*r / J

w1 = w0 - k * 5 *newton*s

pprint("\nω₁ / (1/s):")
tmp = w1
tmp/= 1/s
pprint(N(tmp,prec))

pprint("\ndt / s:")
dt = w1 / k / ( 5 *newton)
tmp = dt
tmp /= s
pprint(tmp)
