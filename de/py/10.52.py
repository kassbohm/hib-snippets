from sympy.physics.units import *
from sympy import *

newton = kg*m/s**2
Nm = newton*m

prec = 5

b = 1 *m
alpha = 6 /s**2
mass = 250 *kg
mu = 20 *kg/m

pprint("\n10.52:")
pprint("\nJzz / (kg m²):")
Jzz = mass*b**2/2 + (2*b)**2*mass
tmp = Jzz
tmp /= kg *m**2
pprint(N(tmp,prec))

pprint("\n10.53:")
pprint("\nJzz / (kg m²):")
m1 = mu*(2*b)
m2 = mu*b
Jzz += m1*(2*b)**2/3 + m2*(2*b)**2
tmp = Jzz
tmp /= kg *m**2
pprint(N(tmp,prec))

pprint("\nJzz wz' / (Nm):")
tmp = Jzz*alpha
tmp /= Nm
pprint(N(tmp,prec))

# 10.52:

# Jzz / (kg m²):
# 1125.0

# 10.53:
#
# Jzz / (kg m²):
# 1258.3
#
# Jzz wz' / (Nm):
# 7550.0
