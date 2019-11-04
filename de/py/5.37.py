from sympy.physics.units import *
from sympy import *

r, R, phi, omega = var("r, R, phi, omega")

cp, sp = cos(phi), sin(phi)

xC, yC = R*cp, R*sp

s = sqrt((r+R)**2 - yC**2)

xD = xC + s

pprint("\nSolution Kai:")
pprint("\nvD / w R:")
tmp = diff(xD, phi)*omega
tmp /= omega*R
tmp = tmp.simplify()
pprint(tmp)

pprint("\nHibbeler's solution is identical:")
pprint("\nvD / w R:")
tmp = - R*R*omega*sin(2*phi)
tmp /= 2*sqrt(R*R*cp*cp + r*r+2*r*R)
tmp -= R*omega*sp
tmp /= omega*R
tmp = tmp.simplify()
pprint(tmp)

# Solution Kai:
#
# vD / w R:
#               √2⋅R⋅sin(2⋅φ)
# - ────────────────────────────────────── - sin(φ)
#        _________________________________
#       ╱  2             2              2
#   2⋅╲╱  R ⋅cos(2⋅φ) + R  + 4⋅R⋅r + 2⋅r
#
# Hibbeler's solution is identical:
#
# vD / w R:
#               √2⋅R⋅sin(2⋅φ)
# - ────────────────────────────────────── - sin(φ)
#        _________________________________
#       ╱  2             2              2
#   2⋅╲╱  R ⋅cos(2⋅φ) + R  + 4⋅R⋅r + 2⋅r
