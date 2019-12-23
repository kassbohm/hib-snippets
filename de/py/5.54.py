from sympy.physics.units import *
from sympy import *

wAB, theta, lAB = var("wAB, theta, lAB")
ct, st = cos(theta), sin(theta)
rAB = Matrix([lAB*ct, lAB*st, 0])
WAB = Matrix([0,0,wAB])

vB = WAB.cross(rAB)
pprint(vB)