from sympy.physics.units import *
from sympy import *

( b, c, d ) = ( S(5)/4 *m, 1 *m, S(7)/4 *m)

# Omega and Omega':
O  = Matrix([0, 0, 5]) / s
Op = Matrix([0, 0, 2]) / (s**2)
w  = Matrix([0, 2, 0]) / s
wp = Matrix([0, 1, 0]) / (s**2)

dAB = Matrix([-b, 0, 0])
dBC = Matrix([0, d, c])

vB = O.cross(dAB)
aB = Op.cross(dAB) + O.cross(O.cross(dAB))

# units:
vu = m/s
au = m/(s**2)

# precision:
prec = 10

pprint("\nvB / (m/s):")
pprint(N(vB/vu,prec))
pprint("\naB / (m/s²):")
pprint(N(aB/au,prec))

vCB = w.cross(dBC)
aCB = wp.cross(dBC) + w.cross(w.cross(dBC))

vC = vB + O.cross(dBC) + vCB
aC = aB + Op.cross(dBC) + O.cross(O.cross(dBC)) + 2*O.cross(vCB) + aCB

pprint("\nvC / (m/s):")
pprint(N(vC/vu,prec))
pprint("\naC / (m/s²):")
pprint(N(aC/au,prec))
