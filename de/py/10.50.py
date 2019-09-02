from sympy.physics.units import *
from sympy import *

prec = 5

(b, c, d) = (0.5 *m , 0.3 *m, 0.4 *m)
mu = 1.5 *kg/m

pprint("\nJxy / (kg m²):")
Jxy = - mu*b*d**2/2
tmp = Jxy
tmp /= kg *m**2
pprint(N(tmp,prec))

pprint("\nJyz / (kg m²):")
Jyz = - mu*(b+c)*c**2/2
tmp = Jyz
tmp /= kg *m**2
pprint(N(tmp,prec))

pprint("\nJyy / (kg m²):")
Jyy = mu/3*(d**3 + c**3)
tmp = Jyy
tmp /= kg *m**2
pprint(N(tmp,prec))

# Jxy / (kg m²):
# -0.060000
#
# Jyz / (kg m²):
# -0.054000
#
# Jyy / (kg m²):
# 0.045500
