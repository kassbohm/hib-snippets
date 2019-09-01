from sympy.physics.units import kg, m, s
from sympy import *

newton = kg*m/s**2
kN = 1000*newton

mass = 150 *kg
b = S(4)/10 *m
F = S(8)/10 *kN

prec = 5

pprint("\n\nI wrt axis parallel to edges through CG / (kg m²):")
b2 = b**2
I = mass/12*Matrix([
    [9*b2,  0,  0],
    [0, b2 + 9*b2, 0],
    [0, 0, b2],
    ])
tmp = I
tmp /= (kg * m**2)
pprint(N(tmp,prec))

pprint("\n\nphi / deg:")
phi = atan2(1,3)
tmp = phi*180/pi
pprint(N(tmp,prec))
C, S = cos(phi), sin(phi)

# Unit vector:
u = Matrix([S, 0, C])

# Mom. of inertia wrt axis along u:
pprint("\n\nIuu / (kg m²):")
Iuu = u.transpose().dot(I.dot(u))
tmp = Iuu
tmp /= (kg * m**2)
pprint(N(tmp,prec))

w = var("omega", positive=True)
w2 = w**2
Ekin = Iuu*w2/2
A = F * 3*b*S * 2*pi

eq = Eq(Ekin, A)
sol = solve(eq, w)
pprint("\n\nw / (1/s):")
tmp = sol[0]/(1/s)
pprint(N(tmp,prec))
