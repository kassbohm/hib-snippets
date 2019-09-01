from sympy.physics.units import *
from sympy import *

prec = 3

newton=kg*m/s**2
kilo = 1000

g = S(981)/100 *m/s**2

GZ, GH = 20 *kilo*newton, 10 *kilo*newton
d, h = 2 *m, 2 *m
R = d / 2

# Mass, Volume, Density:
mZ, mH = GZ/g, GH/g
VZ, VH = pi*R**2*h, 4*pi*R**3/6
rZ, rH = mZ/VZ, mH/VH

# h, mH, mZ, rH, R = var("h, m_H, m_Z, rho_H, R")

# Spherical Coordinates:
r, the, phi = var("r, theta, phi")
st, ct = sin(the), cos(the)
sp, cp = sin(phi), cos(phi)

# Cartesian to Spherical Coords:
(x, y, z) = (r * st * cp, r * st * sp, r * ct)

# dV = dx dy dz = jac * dr dθ dφ
jac = r**2 * st

# Components wrt (x',y',z'):
jxx = y**2 + z**2
jyy = x**2 + z**2
jzz = x**2 + y**2
jxy = - x*y
jxz = - x*z
jyz = - y*z

# Integrating over volume of semisphere:
# upper:   0 < θ < π/2
# lower semisphere: π/2 < θ < π
for x in [jxx, jyy, jzz, jxy, jxz, jyz]:
    x *= jac*rH
    tmp = integrate(x, (phi, 0, 2*pi))
    # upper:
    tmp = integrate(tmp, (the, 0, pi/2))
    # lower:
    # tmp = integrate(tmp, (the, pi/2, pi))
    tmp = integrate(tmp, (r, 0, R))

Ixx = 4*pi*R**5*rH/15
Iyy = 4*pi*R**5*rH/15
Izz = 4*pi*R**5*rH/15
Ixy = 0
Ixz = 0
Iyz = 0

# Components wrt (x'',y'',z''):
d3 = 3*R/8
Ixx -= d3**2 * mH

# Components wrt (x,y,z):
d3 += h/2
Ixx += d3**2 * mH

pprint("\nUnion of both Semispheres:")
Ixx_2H = Ixx*2
Izz_2H = Izz*2
unit = 1000 * kg * m**2

pprint("\nIxx / (10³ kg m²):")
pprint(N(Ixx_2H / unit,prec))
pprint("\nIzz / (10³ kg m²):")
pprint(N(Izz_2H / unit ,prec))

pprint("\nCylinder:")
Ixx_Z = mZ*(3*R*R + h*h)/12
Izz_Z = mZ*R**2/2
pprint("\nIxx / (10³ kg m²):")
pprint(N(Ixx_Z / unit ,prec))
pprint("\nIzz / (10³ kg m²):")
pprint(N(Izz_Z / unit,prec))

pprint("\nTotal:")
pprint("\nIxx / (10³ kg m²):")
Ixx = Ixx_2H + Ixx_Z
tmp = Ixx / unit
pprint(N(tmp,prec))
pprint("\nIzz / (10³ kg m²):")
Izz = Izz_2H + Izz_Z
tmp = Izz / unit
pprint(N(tmp,prec))

pprint("\nTotal mass / (10³ kg):")
mass = 2*mH + mZ
tmp = mass
tmp /= 1000*kg
pprint(N(tmp,prec))

I = Matrix([
[Ixx, Ixy, Ixz],
[Ixy, Ixx, Iyz],
[Ixy, Iyz, Izz],
])
tmp = I / unit
tmp = N(tmp,prec)
pprint(u"\nComponents of I wrt (x,y,z) / (10³ kg m²):")
pprint(tmp)

a = var("alpha")

ca, sa = cos(a), sin(a)

e = Matrix([0, ca, sa])
eT = e.transpose()

# with unicode overline=overbar:
pprint("\nI"+u'z\u0304z\u0304'+" / (10³ kg m²):")
tmp = eT.dot(I.dot(e))
tmp = tmp.subs(a, 45 *pi/180)
tmp /= unit
pprint(N(tmp,prec))
