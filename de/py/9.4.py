l = 1.4 *m

theta = N(25*pi/180)

ct, st = cos(theta), sin(theta)

# Omega and Omega' from (1) and (2):
w  = Matrix([2, 0, 6]) / s
wp = Matrix([1.5, 12, 3]) / (s**2)

r = l*Matrix([0, ct, st])

vA = w.cross(r)
aA = wp.cross(r) + w.cross(w.cross(r))

pprint("\nvA / (m/s):")
tmp = vA / (m/s)
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\naA / (m/s²):")
tmp = aA / (m/(s**2))
tmp = iso_round(tmp,0.1)
pprint(tmp)

# vA / (m/s):
# ⎡-7.6⎤
# ⎢    ⎥
# ⎢-1.2⎥
# ⎢    ⎥
# ⎣2.5 ⎦
#
# aA / (m/s²):
# ⎡10.4 ⎤
# ⎢     ⎥
# ⎢-51.6⎥
# ⎢     ⎥
# ⎣-0.5 ⎦
