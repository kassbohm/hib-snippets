from sympy.physics.units import kg, m, s
from sympy import *

# 10.33 and 10.34 solved here
newton = kg*m/s**2

mass = 15 *kg
r = S(8)/10 *m
phi = pi/4
w1 = 8 /s
bb, cc = 4 *newton*m, S(1)/10 /s

# ---

prec = 3
case = 10.33
case = 10.34

half, fourth = S(1)/2, S(1)/4

t, w2, psi = var("t, w2, psi")

# start and end time:

if case == 10.33:
    M = 2 *newton*m
elif case == 10.34:
    M = bb*exp(cc*t)
else:
    pprint("\nError")
    exit()

# xi = t / s:
xi = var("xi")
# Moment in Nm:
pprint("\n\nM(xi) / Nm:")
M_inNm = M / (newton*m)
Mxi = M_inNm.subs(t, xi*s)
pprint(Mxi)

pprint("\n\nM integriert / (Nm s):")
if case == 10.33:
    Mt = integrate(Mxi,(xi,0,3))
elif case == 10.34:
    Mt = integrate(Mxi,(xi,0,2))

tmp = Mt
pprint(N(tmp,prec))


pprint("\n\nI wrt (x,y,z) / (kg m²):")
I = mass*r*r*Matrix([
    [half,  0,  0],
    [0, fourth, 0],
    [0, 0, fourth],
    ])
tmp = I
tmp /= (kg * m**2)
pprint(N(tmp,3))

# Unit vector along AB:
uAB = Matrix([1, 1, 0])
uAB /= uAB.norm()

# Mom. of inertia wrt horizontal AB-axis:
pprint("\n\nIx'x' / (kg m²):")
Ih = uAB.transpose().dot(I.dot(uAB))
tmp = Ih
tmp /= (kg * m**2)
pprint(N(tmp,3))

pprint("\n\nStart:")
pprint("w / (1/s):")
w0 = 8 /s
tmp = w0
tmp *= s
pprint(N(tmp,prec))

pprint("\n\nH / (Nm s):")
H = Ih*w0
tmp = H
tmp /= newton*m*s
pprint(N(tmp,prec))

pprint("\n\nEnd:")
pprint("H / (Nm s):")
H += Mt*newton*m*s
tmp = H
tmp /= newton*m*s
pprint(N(tmp,prec))

w1 = H/Ih
tmp = w1
tmp *= s
pprint("\n\nw / (1/s):")
pprint(N(tmp,prec))
