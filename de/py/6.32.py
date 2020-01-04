g = S(981)/10 *m/s**2

t =        5 *s
d = S(18)/10 *m
G =        1 *kN

mass = G/g

a = var("a")
# x(t) = 1/2 a*t²
eq = Eq(2*d, a*t**2/2)

a = solve(eq,a)[0]
tmp = a
tmp /= m/s**2
tmp = iso_round(tmp,0.001)
pprint(tmp)