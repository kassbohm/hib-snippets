from sympy.physics.units import *
from sympy import *

m, h, r = var("m, h, r")

xi = 3*h/4

Iy = 3*m/20*(r*r + 4*h*h)

pprint("\nIyS")
IyS = Iy - xi*xi*m
tmp = IyS
tmp = tmp.simplify()
pprint(tmp)
