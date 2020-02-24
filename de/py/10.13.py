a = 0.3 *m
q = 25 *Newton/m
grav = 9.81 *m/s**2

pprint("\nMass of one part / kg:")
G = a*q
mass = G/grav
tmp = mass
tmp /= kg
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nx-position of centroid / m:")
(xS, yS) = (-2*a/3, a/2)
tmp = xS
tmp /= m
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\ny-position of centroid / m:")
tmp = yS
tmp /= m
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nIx'x' / (kg m²):")
Ix1 = mass * yS*yS
Ix2 = mass * a*a / 12
Ix3 = Ix1
Ix = Ix1 + Ix2 + Ix3
tmp = Ix
tmp /= kg*m**2
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nIy'y' / (kg m²):")
Iy1 = mass * a*a / 12 + (a/6)**2*mass
Iy2 = 0 + (a/3)**2*mass
Iy3 = Iy1
Iy = Iy1+Iy2+Iy3
tmp = Iy
tmp /= kg*m**2
tmp = iso_round(tmp,0.001)
pprint(tmp)
